// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
#include <cmath>

#include "Eigen/Eigenvalues"
#include "boost/math/special_functions/erf.hpp"

#include "ndarray/eigen.h"

#include "lsst/log/Log.h"
#include "lsst/afw/table/Catalog.h"
#include "lsst/afw/table/BaseTable.h"
#include "lsst/afw/table/BaseRecord.h"
#include "lsst/meas/modelfit/optimizer.h"
#include "lsst/meas/modelfit/Likelihood.h"
#include "lsst/meas/modelfit/Prior.h"

namespace lsst { namespace meas { namespace modelfit {

// ----------------- OptimizerObjective ---------------------------------------------------------------------

void OptimizerObjective::fillObjectiveValueGrid(
    ndarray::Array<Scalar const,2,1> const & grid,
    ndarray::Array<Scalar,1,1> const & output
) const {
    ndarray::Array<Scalar,1,1> residuals = ndarray::allocate(dataSize);
    for (int i = 0, n = output.getSize<0>(); i < n; ++i) {
        computeResiduals(grid[i], residuals);
        output[i] = 0.5*ndarray::asEigenMatrix(residuals).squaredNorm();
        if (hasPrior()) {
            Scalar prior = computePrior(grid[i]);
            output[i] -= std::log(prior);
            if (std::isnan(output[i])) {
                output[i] = std::numeric_limits<Scalar>::infinity();
            }
        }
    }
}

namespace {

class LikelihoodOptimizerObjective : public OptimizerObjective {
public:

    LikelihoodOptimizerObjective(std::shared_ptr<Likelihood> likelihood, std::shared_ptr<Prior> prior) :
        OptimizerObjective(
            likelihood->getDataDim(), likelihood->getNonlinearDim() + likelihood->getAmplitudeDim()
        ),
        _likelihood(likelihood), _prior(prior),
        _modelMatrix(ndarray::allocate(likelihood->getDataDim(), likelihood->getAmplitudeDim()))
    {}

    void computeResiduals(
        ndarray::Array<Scalar const,1,1> const & parameters,
        ndarray::Array<Scalar,1,1> const & residuals
    ) const override {
        int nlDim = _likelihood->getNonlinearDim();
        int ampDim = _likelihood->getAmplitudeDim();
        _likelihood->computeModelMatrix(_modelMatrix, parameters[ndarray::view(0, nlDim)]);
        ndarray::asEigenMatrix(residuals) = ndarray::asEigenMatrix(_modelMatrix).cast<Scalar>()
            * ndarray::asEigenMatrix(parameters[ndarray::view(nlDim, nlDim+ampDim)]);
        auto likelihoodData = _likelihood->getData();
        ndarray::asEigenMatrix(residuals) -= ndarray::asEigenMatrix(likelihoodData).cast<Scalar>();
    }

    bool hasPrior() const override { return static_cast<bool>(_prior); }

    Scalar computePrior(ndarray::Array<Scalar const,1,1> const & parameters) const override {
        int nlDim = _likelihood->getNonlinearDim();
        int ampDim = _likelihood->getAmplitudeDim();
        return _prior->evaluate(parameters[ndarray::view(0, nlDim)],
                                parameters[ndarray::view(nlDim, nlDim+ampDim)]);
    }

    void differentiatePrior(
        ndarray::Array<Scalar const,1,1> const & parameters,
        ndarray::Array<Scalar,1,1> const & gradient,
        ndarray::Array<Scalar,2,1> const & hessian
    ) const override {
        int nlDim = _likelihood->getNonlinearDim();
        int ampDim = _likelihood->getAmplitudeDim();
        int totDim = nlDim + ampDim;
        _prior->evaluateDerivatives(
            parameters[ndarray::view(0, nlDim)],
            parameters[ndarray::view(nlDim, totDim)],
            gradient[ndarray::view(0, nlDim)],
            gradient[ndarray::view(nlDim, totDim)],
            hessian[ndarray::view(0, nlDim)(0, nlDim)],
            hessian[ndarray::view(nlDim, totDim)(nlDim, totDim)],
            hessian[ndarray::view(0, nlDim)(nlDim, totDim)]
        );
    }

private:
    std::shared_ptr<Likelihood> _likelihood;
    std::shared_ptr<Prior> _prior;
    ndarray::Array<Pixel,2,-1> _modelMatrix;
};

} // anonymous

std::shared_ptr<OptimizerObjective> OptimizerObjective::makeFromLikelihood(
    std::shared_ptr<Likelihood> likelihood,
    std::shared_ptr<Prior> prior
) {
    return std::make_shared<LikelihoodOptimizerObjective>(likelihood, prior);
}

// ----------------- Optimizer::IterationData -----------------------------------------------------------------

Optimizer::IterationData::IterationData(int dataSize, int parameterSize) :
    objectiveValue(0.0), priorValue(0.0),
    parameters(ndarray::allocate(parameterSize)),
    residuals(ndarray::allocate(dataSize))
{}

void Optimizer::IterationData::swap(IterationData & other) {
    std::swap(objectiveValue, other.objectiveValue);
    std::swap(priorValue, other.priorValue);
    parameters.swap(other.parameters);
    residuals.swap(other.residuals);
}

// ----------------- OptimizerHistoryRecorder ---------------------------------------------------------------

OptimizerHistoryRecorder::OptimizerHistoryRecorder(
    afw::table::Schema & schema,
    std::shared_ptr<Model> model,
    bool doSaveDerivatives
) :
    outer(
        schema.addField(afw::table::Field<int>("outer", "current outer iteration count"), true)
    ),
    inner(
        schema.addField(afw::table::Field<int>("inner", "current inner iteration count"), true)
    ),
    state(
        schema.addField(
            afw::table::Field<int>(
                "state", "state bitflags after this step; see Optimizer::StateFlags"
            ),
            true
        )
    ),
    objective(
        schema.addField(
            afw::table::Field<Scalar>(
                "objective", "value of objective function (-ln P) at parameters"
            ),
            true
        )
    ),
    prior(
        schema.addField(afw::table::Field<Scalar>("prior", "prior probability at parameters"), true)
    ),
    trust(
        schema.addField(afw::table::Field<Scalar>("trust", "size of trust region after this step"), true)
    ),
    parameters(
        schema.addField(
            afw::table::Field<afw::table::Array<Scalar> >(
                "parameters",
                "parameter vector",
                model->getNonlinearDim() + model->getAmplitudeDim()
            ),
            true
        )
    )
{
    if (doSaveDerivatives) {
        int const n = model->getNonlinearDim() + model->getAmplitudeDim();
        derivatives = schema.addField(
            afw::table::Field<afw::table::Array<Scalar> >(
                "derivatives",
                "objective function derivatives; use unpackDerivatives() to unpack",
                n + n*(n+1)/2
            ),
            true
        );
    }
}

OptimizerHistoryRecorder::OptimizerHistoryRecorder(afw::table::Schema const & schema) :
    outer(schema["outer"]),
    inner(schema["inner"]),
    state(schema["state"]),
    objective(schema["objective"]),
    prior(schema["prior"]),
    trust(schema["trust"]),
    parameters(schema["parameters"])
{
    try {
        derivatives = schema["derivatives"];
    } catch (pex::exceptions::NotFoundError &) {}
}

void OptimizerHistoryRecorder::apply(
    int outerIterCount,
    int innerIterCount,
    afw::table::BaseCatalog & history,
    Optimizer const & optimizer
) const {
    std::shared_ptr<afw::table::BaseRecord> record = history.addNew();
    record->set(outer, outerIterCount);
    record->set(inner, innerIterCount);
    record->set(state, optimizer.getState());
    record->set(trust, optimizer._trustRadius);
    Optimizer::IterationData const * data;
    if (!(optimizer.getState() & Optimizer::STATUS_STEP_REJECTED)) {
        data = &optimizer._current;
        if (derivatives.isValid()) {
            int const n = parameters.getSize();
            ndarray::Array<Scalar,1,1> packed = (*record)[derivatives];
            for (int i = 0, k = n; i < n; ++i) {
                packed[i] = optimizer._gradient[i];
                for (int j = 0; j <= i; ++j, ++k) {
                    packed[k] = optimizer._hessian(i, j);
                }
            }
        }
    } else {
        data = &optimizer._next;
    }
    record->set(parameters, data->parameters);
    record->set(objective, data->objectiveValue);
    record->set(prior, data->priorValue);
}

void OptimizerHistoryRecorder::unpackDerivatives(
    ndarray::Array<Scalar const,1,1> const & packed,
    Vector & gradient,
    Matrix & hessian
) const {
    int const n = parameters.getSize();
    for (int i = 0, k = n; i < n; ++i) {
        gradient[i] = packed[i];
        for (int j = 0; j <= i; ++j, ++k) {
            hessian(i, j) = hessian(j, i) = packed[k];
        }
    }
}

void OptimizerHistoryRecorder::unpackDerivatives(
    afw::table::BaseRecord const & record,
    ndarray::Array<Scalar,1,1> const & gradient,
    ndarray::Array<Scalar,2,2> const & hessian
) const {
    if (!derivatives.isValid()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicError,
            "HistoryRecorder was not configured to save derivatives"
        );
    }
    return unpackDerivatives(record[derivatives], gradient, hessian);
}

void OptimizerHistoryRecorder::unpackDerivatives(
    ndarray::Array<Scalar const,1,1> const & packed,
    ndarray::Array<Scalar,1,1> const & gradient,
    ndarray::Array<Scalar,2,2> const & hessian
) const {
    int const n = parameters.getSize();
    for (int i = 0, k = n; i < n; ++i) {
        gradient[i] = packed[i];
        for (int j = 0; j <= i; ++j, ++k) {
            hessian[i][j] = hessian[j][i] = packed[k];
        }
    }
}

void OptimizerHistoryRecorder::unpackDerivatives(
    afw::table::BaseRecord const & record,
    Vector & gradient,
    Matrix & hessian
) const {
    if (!derivatives.isValid()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicError,
            "HistoryRecorder was not configured to save derivatives"
        );
    }
    return unpackDerivatives(record[derivatives], gradient, hessian);
}

void OptimizerHistoryRecorder::fillObjectiveModelGrid(
    afw::table::BaseRecord const & record,
    ndarray::Array<Scalar const,2,1> const & grid,
    ndarray::Array<Scalar,1,1> const & output
) const {
    Scalar q = record.get(objective);
    Vector gradient(parameters.getSize());
    Matrix hessian(parameters.getSize(), parameters.getSize());
    Vector s(parameters.getSize());
    // currentNdArray must be a local variable because it owns the data in `current`
    auto currentNdArray = record.get(parameters);
    Vector current = ndarray::asEigenMatrix(currentNdArray);
    unpackDerivatives(record, gradient, hessian);
    for (int i = 0, n = output.getSize<0>(); i < n; ++i) {
        s = ndarray::asEigenMatrix(grid[i]) - current;
        output[i] = q + s.dot(gradient + 0.5*hessian*s);
    }
}

// ----------------- Optimizer ------------------------------------------------------------------------------

Optimizer::Optimizer(
    std::shared_ptr<Objective const> objective,
    ndarray::Array<Scalar const,1,1> const & parameters,
    Control const & ctrl
) :
    _state(0x0),
    _objective(objective),
    _ctrl(ctrl),
    _trustRadius(ctrl.trustRegionInitialSize),
    _current(objective->dataSize, objective->parameterSize),
    _next(objective->dataSize, objective->parameterSize),
    _step(ndarray::allocate(objective->parameterSize)),
    _gradient(ndarray::allocate(objective->parameterSize)),
    _hessian(ndarray::allocate(objective->parameterSize, objective->parameterSize)),
    _residualDerivative(ndarray::allocate(objective->dataSize, objective->parameterSize)),
    _sr1b(objective->parameterSize, objective->parameterSize),
    _sr1v(objective->parameterSize),
    _sr1jtr(objective->parameterSize)
{
    LOG_LOGGER trace3Logger = LOG_GET("TRACE3.meas.modelfit.optimizer.Optimizer");
    if (parameters.getSize<0>() != static_cast<std::size_t>(_objective->parameterSize)) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthError,
            (boost::format("Parameter vector size (%d) does not match objective (%d)")
             % parameters.getSize<0>() % _objective->parameterSize).str()
        );
    }
    _current.parameters.deep() = parameters;
    _next.parameters.deep() = parameters;
    _objective->computeResiduals(_current.parameters, _current.residuals);
    _current.objectiveValue = 0.5*ndarray::asEigenMatrix(_current.residuals).squaredNorm();
    if (_objective->hasPrior()) {
        _current.priorValue = _objective->computePrior(_current.parameters);
        _current.objectiveValue -= std::log(_current.priorValue);
    }
    LOGL_DEBUG(trace3Logger, "Initial objective value is %g", _current.objectiveValue);
    _sr1b.setZero();
    _computeDerivatives();
    ndarray::asEigenMatrix(_hessian) = ndarray::asEigenMatrix(_hessian).selfadjointView<Eigen::Lower>();
}

void Optimizer::_computeDerivatives() {
    auto resDer = ndarray::asEigenMatrix(_residualDerivative);
    resDer.setZero();
    _next.parameters.deep() = _current.parameters;
    if (!_objective->differentiateResiduals(_current.parameters, _residualDerivative)) {
        for (int n = 0; n < _objective->parameterSize; ++n) {
            double numDiffStep = _ctrl.numDiffRelStep * _next.parameters[n]
                + _ctrl.numDiffTrustRadiusStep * _trustRadius
                + _ctrl.numDiffAbsStep;
            _next.parameters[n] += numDiffStep;
            _objective->computeResiduals(_next.parameters, _next.residuals);
            resDer.col(n) =
                    (ndarray::asEigenMatrix(_next.residuals) - ndarray::asEigenMatrix(_current.residuals)) /
                    numDiffStep;
            _next.parameters[n] = _current.parameters[n];
        }
    }
    _gradient.deep() = 0.0;
    _hessian.deep() = 0.0;
    if (_objective->hasPrior()) {
        _objective->differentiatePrior(_current.parameters, _gradient, _hessian);
        // objective evaluates P(x); we want -ln P(x) and associated derivatives
        ndarray::asEigenMatrix(_gradient) /= -_current.priorValue;
        ndarray::asEigenMatrix(_hessian) /= -_current.priorValue;
        ndarray::asEigenMatrix(_hessian).selfadjointView<Eigen::Lower>().rankUpdate(
                ndarray::asEigenMatrix(_gradient), 1.0);
    }
    if (!_ctrl.noSR1Term) {
        _sr1jtr = resDer.adjoint() * ndarray::asEigenMatrix(_current.residuals);
        ndarray::asEigenMatrix(_gradient) += _sr1jtr;
    } else {
        ndarray::asEigenMatrix(_gradient) += resDer.adjoint() * ndarray::asEigenMatrix(_current.residuals);
    }
    ndarray::asEigenMatrix(_hessian).selfadjointView<Eigen::Lower>().rankUpdate(resDer.adjoint(), 1.0);
}

void Optimizer::removeSR1Term() {
   ndarray::asEigenMatrix(_hessian) -= _sr1b;
}

bool Optimizer::_stepImpl(
    int outerIterCount,
    HistoryRecorder const * recorder,
    afw::table::BaseCatalog * history
) {
    LOG_LOGGER trace5Logger = LOG_GET("TRACE5.meas.modelfit.optimizer.Optimizer");
    LOG_LOGGER trace3Logger = LOG_GET("TRACE3.meas.modelfit.optimizer.Optimizer");
    _state &= ~int(STATUS);
    if (ndarray::asEigenMatrix(_gradient).lpNorm<Eigen::Infinity>() <= _ctrl.gradientThreshold) {
        LOGL_DEBUG(trace3Logger, "max(gradient)=%g below threshold; declaring convergence",
                   ndarray::asEigenMatrix(_gradient).lpNorm<Eigen::Infinity>());
        _state |= CONVERGED_GRADZERO;
        return false;
    }
    for (int innerIterCount = 0; innerIterCount < _ctrl.maxInnerIterations; ++innerIterCount) {
        LOGL_DEBUG(trace5Logger, "Starting inner iteration %d", innerIterCount);
        _state &= ~int(STATUS);
        _next.objectiveValue = 0.0;
        _next.priorValue = 1.0;
        solveTrustRegion(
            _step, _hessian, _gradient, _trustRadius, _ctrl.trustRegionSolverTolerance
        );
        ndarray::asEigenMatrix(_next.parameters) =
                ndarray::asEigenMatrix(_current.parameters) + ndarray::asEigenMatrix(_step);
        double stepLength = ndarray::asEigenMatrix(_step).norm();
        if (std::isnan(stepLength)) {
            LOGL_DEBUG(trace3Logger, "NaN encountered in step length");
            _state |= FAILED_NAN;
            return false;
        }
        LOGL_DEBUG(trace5Logger, "Step has length %g", stepLength);
        if (_objective->hasPrior()) {
            _next.priorValue = _objective->computePrior(_next.parameters);
            _next.objectiveValue = -std::log(_next.priorValue);
            if (_next.priorValue <= 0.0 || std::isnan(_next.objectiveValue)) {
                _next.objectiveValue = std::numeric_limits<Scalar>::infinity();
                LOGL_DEBUG(trace5Logger, "Rejecting step due to zero prior");
                if (stepLength < _trustRadius) {
                    // Because the step failed due to the prior, we could add an API to the objective
                    // for projecting the step back to a feasible value, and move directly there.
                    // It remains to be seen if that's needed.
                    LOGL_DEBUG(trace5Logger, "Unconstrained step failed; setting trust radius to step length %g",
                               stepLength);
                    _trustRadius = stepLength;
                }
                _trustRadius *= _ctrl.trustRegionShrinkFactor;
                LOGL_DEBUG(trace5Logger, "Decreasing trust radius to %g", _trustRadius);
                _state |= STATUS_STEP_REJECTED | STATUS_TR_DECREASED;
                if (_trustRadius <= _ctrl.minTrustRadiusThreshold) {
                    LOGL_DEBUG(trace3Logger, "Trust radius %g has dropped below threshold %g; declaring convergence",
                                 _trustRadius, _ctrl.minTrustRadiusThreshold);
                    _state |= CONVERGED_TR_SMALL;
                    return false;
                }
                if (recorder) recorder->apply(outerIterCount, innerIterCount, *history, *this);
                continue;
            }
        }
        _objective->computeResiduals(_next.parameters, _next.residuals);
        _next.objectiveValue += 0.5*ndarray::asEigenMatrix(_next.residuals).squaredNorm();
        double actualChange = _next.objectiveValue - _current.objectiveValue;
        double predictedChange = ndarray::asEigenMatrix(_step).dot(ndarray::asEigenMatrix(_gradient) +
                                                                   0.5 * ndarray::asEigenMatrix(_hessian) *
                                                                           ndarray::asEigenMatrix(_step));
        double rho = actualChange / predictedChange;
        if (std::isnan(rho)) {
            LOGL_DEBUG(trace5Logger, "NaN encountered in rho");
            _state |= FAILED_NAN;
            return false;
        }
        LOGL_DEBUG(trace5Logger, "Reduction ratio rho=%g; actual=%g, predicted=%g", rho, actualChange, predictedChange);
        if (rho > _ctrl.stepAcceptThreshold && actualChange < 0.0) {
            LOGL_DEBUG(trace5Logger, "Step accepted; new objective=%g, old was %g", _next.objectiveValue,
                       _current.objectiveValue);
            _state |= STATUS_STEP_ACCEPTED;
            _current.swap(_next);
            if (!_ctrl.noSR1Term) {
                _sr1v = -_sr1jtr;
            }
            _computeDerivatives();
            if (!_ctrl.noSR1Term) {
                _sr1v += _sr1jtr;
                double vs = _sr1v.dot(ndarray::asEigenMatrix(_step));
                if (vs >= (_ctrl.skipSR1UpdateThreshold * _sr1v.norm() * stepLength + 1.0)) {
                    _sr1b.selfadjointView<Eigen::Lower>().rankUpdate(_sr1v, 1.0 / vs);
                }
                ndarray::asEigenMatrix(_hessian) += _sr1b;
            }
            ndarray::asEigenMatrix(_hessian) =
                    ndarray::asEigenMatrix(_hessian).selfadjointView<Eigen::Lower>();
            if (rho > _ctrl.trustRegionGrowReductionRatio &&
                (stepLength / _trustRadius) > _ctrl.trustRegionGrowStepFraction) {
                _state |= STATUS_TR_INCREASED;
                _trustRadius *= _ctrl.trustRegionGrowFactor;
                LOGL_DEBUG(trace5Logger, "Increasing trust radius to %g", _trustRadius);
            } else if (rho < _ctrl.trustRegionShrinkReductionRatio) {
                // even though the step was accepted, our quadratic model
                // of the objective function wasn't very accurate, so we
                // decrease the trust region anyway
                _state |= STATUS_TR_DECREASED;
                _trustRadius *= _ctrl.trustRegionShrinkFactor;
                LOGL_DEBUG(trace5Logger, "Decreasing trust radius to %g", _trustRadius);
            } else {
                LOGL_DEBUG(trace5Logger, "Leaving trust radius unchanged at %g", _trustRadius);
                _state |= STATUS_TR_UNCHANGED;
            }
            if (recorder) recorder->apply(outerIterCount, innerIterCount, *history, *this);
            return true;
        }
        _state |= STATUS_STEP_REJECTED;
        LOGL_DEBUG(trace5Logger, "Step rejected; test objective was %g, current is %g", _next.objectiveValue,
                   _current.objectiveValue);
        if (stepLength < _trustRadius) {
            LOGL_DEBUG(trace5Logger, "Unconstrained step failed; setting trust radius to step length %g",
                       stepLength);
            _trustRadius = stepLength;
        }
        // we always decrease the trust radius if the step is rejected - otherwise we'll just
        // produce the same step again
        _state |= STATUS_TR_DECREASED;
        _trustRadius *= _ctrl.trustRegionShrinkFactor;
        LOGL_DEBUG(trace5Logger, "Decreasing trust radius to %g", _trustRadius);
        if (_trustRadius <= _ctrl.minTrustRadiusThreshold) {
            _state |= CONVERGED_TR_SMALL;
            LOGL_DEBUG(trace3Logger, "Trust radius %g has dropped below threshold %g; declaring convergence",
                       _trustRadius, _ctrl.minTrustRadiusThreshold);
            return false;
        }
        if (recorder) recorder->apply(outerIterCount, innerIterCount, *history, *this);
    }
    LOGL_DEBUG(trace3Logger, "Max inner iteration number exceeded");
    _state |= FAILED_MAX_INNER_ITERATIONS;
    return false;
}

int Optimizer::_runImpl(HistoryRecorder const * recorder, afw::table::BaseCatalog * history) {
    LOG_LOGGER trace5Logger = LOG_GET("TRACE5.meas.modelfit.optimizer.Optimizer");
    LOG_LOGGER trace3Logger = LOG_GET("TRACE3.meas.modelfit.optimizer.Optimizer");
    if (recorder) recorder->apply(-1, -1, *history, *this);
    int outerIterCount = 0;
    try {
        for (; outerIterCount < _ctrl.maxOuterIterations; ++outerIterCount) {
            LOGL_DEBUG(trace5Logger, "Starting outer iteration %d", outerIterCount);
            if (!_stepImpl(outerIterCount, recorder, history)) return outerIterCount;
        }
        _state |= FAILED_MAX_OUTER_ITERATIONS;
        LOGL_DEBUG(trace3Logger, "Max outer iteration number exceeded");
    } catch (...) {
        _state |= FAILED_EXCEPTION;
    }
    return outerIterCount;
}


// ----------------- Trust Region solver --------------------------------------------------------------------

void solveTrustRegion(
    ndarray::Array<Scalar,1,1> const & x,
    ndarray::Array<Scalar const,2,1> const & F,
    ndarray::Array<Scalar const,1,1> const & g,
    double r, double tolerance
) {
    static double const ROOT_EPS = std::sqrt(std::numeric_limits<double>::epsilon());
    static int const ITER_MAX = 10;
    LOG_LOGGER trace5Logger = LOG_GET("TRACE5.meas.modelfit.optimizer.Optimizer");
    LOG_LOGGER trace3Logger = LOG_GET("TRACE3.meas.modelfit.optimizer.Optimizer");
    double const r2 = r*r;
    double const r2min = r2 * (1.0 - tolerance) * (1.0 - tolerance);
    double const r2max = r2 * (1.0 + tolerance) * (1.0 + tolerance);
    int const d = g.getSize<0>();
    Eigen::SelfAdjointEigenSolver<Matrix> eigh(ndarray::asEigenMatrix(F));
    double const threshold = ROOT_EPS * eigh.eigenvalues()[d - 1];
    Vector qtg = eigh.eigenvectors().adjoint() * ndarray::asEigenMatrix(g);
    Vector tmp(d);
    double mu = 0.0;
    double xsn = 0.0;
    if (eigh.eigenvalues()[0] >= threshold) {
        LOGL_DEBUG(trace5Logger, "Starting with full-rank matrix");
        tmp = (eigh.eigenvalues().array().inverse() * qtg.array()).matrix();
        ndarray::asEigenMatrix(x) = -eigh.eigenvectors() * tmp;
        xsn = ndarray::asEigenMatrix(x).squaredNorm();
        if (xsn <= r2max) {
            LOGL_DEBUG(trace5Logger, "Ending with unconstrained solution");
            // unconstrained solution is within the constraint; no more work to do
            return;
        }
    } else {
        mu = -eigh.eigenvalues()[0] + 2.0*ROOT_EPS*eigh.eigenvalues()[d - 1];
        tmp = ((eigh.eigenvalues().array() + mu).inverse() * qtg.array()).matrix();
        int n = 0;
        while (eigh.eigenvalues()[++n] < threshold);
        LOGL_DEBUG(trace5Logger, "Starting with %d zero eigenvalue(s) (of %d)", n, d);
        if ((qtg.head(n).array() < ROOT_EPS * ndarray::asEigenMatrix(g).lpNorm<Eigen::Infinity>()).all()) {
            ndarray::asEigenMatrix(x) = -eigh.eigenvectors().rightCols(n) * tmp.tail(n);
            xsn = ndarray::asEigenMatrix(x).squaredNorm();
            if (xsn < r2min) {
                // Nocedal and Wright's "Hard Case", which is actually
                // easier: Q_1^T g is zero (where the columns of Q_1
                // are the eigenvectors that correspond to the
                // smallest eigenvalue \lambda_1), so \mu = -\lambda_1
                // and we can add a multiple of any column of Q_1 to x
                // to get ||x|| == r.  If ||x|| > r, we can find the
                // solution with the usual iteration by increasing \mu.
                double tau = std::sqrt(r*r - ndarray::asEigenMatrix(x).squaredNorm());
                ndarray::asEigenMatrix(x) += tau * eigh.eigenvectors().col(0);
                LOGL_DEBUG(trace5Logger, "Ending; Q_1^T g == 0, and ||x|| < r");
                return;
            }
            LOGL_DEBUG(trace5Logger, "Continuing; Q_1^T g == 0, but ||x|| > r");
        } else {
            ndarray::asEigenMatrix(x) = -eigh.eigenvectors() * tmp;
            xsn = ndarray::asEigenMatrix(x).squaredNorm();
            LOGL_DEBUG(trace5Logger, "Continuing; Q_1^T g != 0, ||x||=%f");
        }
    }
    int nIter = 0;
    while ((xsn < r2min || xsn > r2max) && ++nIter < ITER_MAX) {
        LOGL_DEBUG(trace5Logger, "Iterating at mu=%f, ||x||=%f, r=%f", mu, std::sqrt(xsn), r);
        mu += xsn*(std::sqrt(xsn) / r - 1.0)
            / (qtg.array().square() / (eigh.eigenvalues().array() + mu).cube()).sum();
        tmp = ((eigh.eigenvalues().array() + mu).inverse() * qtg.array()).matrix();
        ndarray::asEigenMatrix(x) = -eigh.eigenvectors() * tmp;
        xsn = ndarray::asEigenMatrix(x).squaredNorm();
    }
    LOGL_DEBUG(trace5Logger, "Ending at mu=%f, ||x||=%f, r=%f", mu, std::sqrt(xsn), r);
    return;
}

}}} // namespace lsst::meas::modelfit
