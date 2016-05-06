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

#define LSST_MAX_DEBUG 10
#include "lsst/pex/logging/Debug.h"
#include "lsst/afw/table/Catalog.h"
#include "lsst/afw/table/BaseTable.h"
#include "lsst/afw/table/BaseRecord.h"
#include "lsst/meas/modelfit/optimizer.h"
#include "lsst/meas/modelfit/Likelihood.h"
#include "lsst/meas/modelfit/Prior.h"
#include "lsst/meas/modelfit/ModelFitRecord.h"

namespace lsst { namespace meas { namespace modelfit {

// ----------------- OptimizerInterpreter -------------------------------------------------------------------

namespace {

Model::NameVector concatenateNameVectors(Model::NameVector const & a, Model::NameVector const & b) {
    Model::NameVector r;
    r.reserve(a.size() + b.size());
    r.insert(r.end(), a.begin(), a.end());
    r.insert(r.end(), b.begin(), b.end());
    return r;
}

ndarray::Array<Scalar,1,1> computeGaussianQuantile(
    Mixture const & pdf, int dim, ndarray::Array<Scalar const,1,1> const & fractions
) {
    MixtureComponent const & component = *pdf.begin();
    // TODO: this could be more efficient if we made changes to MixtureComponent;
    // we don't need to get the full vector+matrix just to get an element of each.
    Scalar mu = component.getMu()[dim];
    Scalar sigma = component.getSigma()(dim, dim);
    ndarray::Array<Scalar,1,1> output = ndarray::allocate(fractions.getSize<0>());
    for (int i = 0, n = fractions.getSize<0>(); i < n; ++i) {
        output[i] = mu + M_SQRT2*sigma*boost::math::erf_inv(2.0*fractions[i] - 1.0);
    }
    return output;
}

} // anonymous

OptimizerInterpreter::OptimizerInterpreter(PTR(Model) model, PTR(Prior) prior) :
    Interpreter(concatenateNameVectors(model->getNonlinearNames(), model->getAmplitudeNames()), model, prior)
{}

void OptimizerInterpreter::attachPdf(ModelFitRecord & record, Optimizer const & optimizer) const {
    static Scalar const THRESHOLD = 1E-8;
    ndarray::Array<Scalar const,2,2> hessian = optimizer.getHessian();
    Eigen::SelfAdjointEigenSolver<Matrix> eig(hessian.asEigen().adjoint());
    Vector s = eig.eigenvalues();
    // fudge eigenvalues to ensure the matrix is invertible; this transforms fully degenerate directions
    // in parameter space into nearly-degenerate directions.
    for (int i = 0, n = s.size()-1; i < n; ++i) {
        if (s[i] < THRESHOLD*s[n]) {
            s[i] = THRESHOLD*s[n];
        }
    }
    Matrix sigma = eig.eigenvectors() * s.array().inverse().matrix().asDiagonal()
        * eig.eigenvectors().adjoint();
    Mixture::ComponentList components;
    components.push_back(
        Mixture::Component(1.0, optimizer.getParameters().asEigen(), sigma)
    );
    record.setPdf(std::make_shared<Mixture>(s.size(), boost::ref(components)));
}

ndarray::Array<Scalar,1,1> OptimizerInterpreter::computeParameterQuantiles(
    ModelFitRecord const & record,
    ndarray::Array<Scalar const,1,1> const & fractions,
    int index
) const {
    // TODO: bounds checking
    return computeGaussianQuantile(*record.getPdf(), index, fractions);
}

ndarray::Array<Scalar,1,1> OptimizerInterpreter::computeNonlinearQuantiles(
    ModelFitRecord const & record,
    ndarray::Array<Scalar const,1,1> const & fractions,
    int index
) const {
    // TODO: bounds checking
    return computeGaussianQuantile(*record.getPdf(), index, fractions);
}

ndarray::Array<Scalar,1,1> OptimizerInterpreter::computeAmplitudeQuantiles(
    ModelFitRecord const & record,
    ndarray::Array<Scalar const,1,1> const & fractions,
    int index
) const {
    // TODO: bounds checking
    return computeGaussianQuantile(*record.getPdf(), index + getModel()->getNonlinearDim(), fractions);
}

ndarray::Array<Scalar,1,1> OptimizerInterpreter::computeParameterMean(ModelFitRecord const & record) const {
    ndarray::Array<Scalar,1,1> output = ndarray::allocate(getParameterDim());
    output.asEigen() = record.getPdf()->begin()->getMu();
    return output;
}

ndarray::Array<Scalar,1,1> OptimizerInterpreter::computeNonlinearMean(ModelFitRecord const & record) const {
    ndarray::Array<Scalar,1,1> output = ndarray::allocate(getNonlinearDim());
    output.asEigen() = record.getPdf()->begin()->getMu().head(getNonlinearDim());
    return output;
}

ndarray::Array<Scalar,1,1> OptimizerInterpreter::computeAmplitudeMean(ModelFitRecord const & record) const {
    ndarray::Array<Scalar,1,1> output = ndarray::allocate(getAmplitudeDim());
    output.asEigen() = record.getPdf()->begin()->getMu().tail(getAmplitudeDim());
    return output;
}

ndarray::Array<Scalar,2,2> OptimizerInterpreter::computeParameterCovariance(
    ModelFitRecord const & record,
    ndarray::Array<Scalar const,1,1> const & mean
) const {
    ndarray::Array<Scalar,2,2> output = ndarray::allocate(getParameterDim(), getParameterDim());
    output.asEigen() = record.getPdf()->begin()->getSigma().adjoint();
    return output;
}

ndarray::Array<Scalar,2,2> OptimizerInterpreter::computeNonlinearCovariance(
    ModelFitRecord const & record,
    ndarray::Array<Scalar const,1,1> const & mean
) const {
    int const n = getNonlinearDim();
    ndarray::Array<Scalar,2,2> output = ndarray::allocate(n, n);
    output.asEigen() = record.getPdf()->begin()->getSigma().topLeftCorner(n, n).adjoint();
    return output;
}

ndarray::Array<Scalar,2,2> OptimizerInterpreter::computeAmplitudeCovariance(
     ModelFitRecord const & record,
     ndarray::Array<Scalar const,1,1> const & mean
) const {
    int const n = getAmplitudeDim();
    ndarray::Array<Scalar,2,2> output = ndarray::allocate(n, n);
    output.asEigen() = record.getPdf()->begin()->getSigma().bottomRightCorner(n, n).adjoint();
    return output;
}

void OptimizerInterpreter::_packParameters(
    ndarray::Array<Scalar const,1,1> const & nonlinear,
    ndarray::Array<Scalar const,1,1> const & amplitudes,
    ndarray::Array<Scalar,1,1> const & parameters
) const {
    parameters[ndarray::view(0, getNonlinearDim())] = nonlinear;
    parameters[ndarray::view(getNonlinearDim(), getParameterDim())] = amplitudes;
}

void OptimizerInterpreter::_unpackNonlinear(
    ndarray::Array<Scalar const,1,1> const & parameters,
    ndarray::Array<Scalar,1,1> const & nonlinear
) const {
    nonlinear.deep() = parameters[ndarray::view(0, getNonlinearDim())];
}

// ----------------- OptimizerObjective ---------------------------------------------------------------------

void OptimizerObjective::fillObjectiveValueGrid(
    ndarray::Array<Scalar const,2,1> const & grid,
    ndarray::Array<Scalar,1,1> const & output
) const {
    ndarray::Array<Scalar,1,1> residuals = ndarray::allocate(dataSize);
    for (int i = 0, n = output.getSize<0>(); i < n; ++i) {
        computeResiduals(grid[i], residuals);
        output[i] = 0.5*residuals.asEigen().squaredNorm();
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

    LikelihoodOptimizerObjective(PTR(Likelihood) likelihood, PTR(Prior) prior) :
        OptimizerObjective(
            likelihood->getDataDim(), likelihood->getNonlinearDim() + likelihood->getAmplitudeDim()
        ),
        _likelihood(likelihood), _prior(prior),
        _modelMatrix(ndarray::allocate(likelihood->getDataDim(), likelihood->getAmplitudeDim()))
    {}

    virtual void computeResiduals(
        ndarray::Array<Scalar const,1,1> const & parameters,
        ndarray::Array<Scalar,1,1> const & residuals
    ) const {
        int nlDim = _likelihood->getNonlinearDim();
        int ampDim = _likelihood->getAmplitudeDim();
        _likelihood->computeModelMatrix(_modelMatrix, parameters[ndarray::view(0, nlDim)]);
        residuals.asEigen() = _modelMatrix.asEigen().cast<Scalar>()
            * parameters[ndarray::view(nlDim, nlDim+ampDim)].asEigen();
        residuals.asEigen() -= _likelihood->getData().asEigen().cast<Scalar>();
    }

    virtual bool hasPrior() const { return static_cast<bool>(_prior); }

    virtual Scalar computePrior(ndarray::Array<Scalar const,1,1> const & parameters) const {
        int nlDim = _likelihood->getNonlinearDim();
        int ampDim = _likelihood->getAmplitudeDim();
        return _prior->evaluate(parameters[ndarray::view(0, nlDim)],
                                parameters[ndarray::view(nlDim, nlDim+ampDim)]);
    }

    virtual void differentiatePrior(
        ndarray::Array<Scalar const,1,1> const & parameters,
        ndarray::Array<Scalar,1,1> const & gradient,
        ndarray::Array<Scalar,2,1> const & hessian
    ) const {
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
    PTR(Likelihood) _likelihood;
    PTR(Prior) _prior;
    ndarray::Array<Pixel,2,-1> _modelMatrix;
};

} // anonymous

PTR(OptimizerObjective) OptimizerObjective::makeFromLikelihood(
    PTR(Likelihood) likelihood,
    PTR(Prior) prior
) {
    return std::make_shared<LikelihoodOptimizerObjective>(likelihood, prior);
}

Scalar OptimizerObjective::computePrior(ndarray::Array<Scalar const,1,1> const & parameters) const {
    return 1.0;
}

void OptimizerObjective::differentiatePrior(
    ndarray::Array<Scalar const,1,1> const & parameters,
    ndarray::Array<Scalar,1,1> const & gradient,
    ndarray::Array<Scalar,2,1> const & hessian
) const {
    gradient.deep() = 0.0;
    hessian.deep() = 0.0;
}

// ----------------- OptimizerIterationData -----------------------------------------------------------------

OptimizerIterationData::OptimizerIterationData(int dataSize, int parameterSize) :
    objectiveValue(0.0), priorValue(0.0),
    parameters(ndarray::allocate(parameterSize)),
    residuals(ndarray::allocate(dataSize))
{}

void OptimizerIterationData::swap(OptimizerIterationData & other) {
    std::swap(objectiveValue, other.objectiveValue);
    std::swap(priorValue, other.priorValue);
    parameters.swap(other.parameters);
    residuals.swap(other.residuals);
}

// ----------------- OptimizerHistoryRecorder ---------------------------------------------------------------

OptimizerHistoryRecorder::OptimizerHistoryRecorder(
    afw::table::Schema & schema,
    PTR(Model) model,
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
    PTR(afw::table::BaseRecord) record = history.addNew();
    record->set(outer, outerIterCount);
    record->set(inner, innerIterCount);
    record->set(state, optimizer.getState());
    record->set(trust, optimizer._trustRadius);
    OptimizerIterationData const * data;
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
    Vector current = record.get(parameters).asEigen();
    unpackDerivatives(record, gradient, hessian);
    for (int i = 0, n = output.getSize<0>(); i < n; ++i) {
        s = grid[i].asEigen() - current;
        output[i] = q + s.dot(gradient + 0.5*hessian*s);
    }
}

// ----------------- Optimizer ------------------------------------------------------------------------------

Optimizer::Optimizer(
    PTR(Objective const) objective,
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
    _jacobian(objective->dataSize, objective->parameterSize),
    _sr1b(objective->parameterSize, objective->parameterSize),
    _sr1v(objective->parameterSize),
    _sr1jtr(objective->parameterSize)
{
    pex::logging::Debug log("meas.modelfit.optimizer.Optimizer");
    if (parameters.getSize<0>() != _objective->parameterSize) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthError,
            (boost::format("Parameter vector size (%d) does not match objective (%d)")
             % parameters.getSize<0>() % _objective->parameterSize).str()
        );
    }
    _current.parameters.deep() = parameters;
    _next.parameters.deep() = parameters;
    _objective->computeResiduals(_current.parameters, _current.residuals);
    _current.objectiveValue = 0.5*_current.residuals.asEigen().squaredNorm();
    if (_objective->hasPrior()) {
        _current.priorValue = _objective->computePrior(_current.parameters);
        _current.objectiveValue -= std::log(_current.priorValue);
    }
    log.debug<7>("Initial objective value is %g", _current.objectiveValue);
    _sr1b.setZero();
    _computeDerivatives();
    _hessian.asEigen() = _hessian.asEigen().selfadjointView<Eigen::Lower>();
}

void Optimizer::_computeDerivatives() {
    _jacobian.setZero();
    _next.parameters.deep() = _current.parameters;
    for (int n = 0; n < _objective->parameterSize; ++n) {
        double numDiffStep = _ctrl.numDiffRelStep * _next.parameters[n]
            + _ctrl.numDiffTrustRadiusStep * _trustRadius
            + _ctrl.numDiffAbsStep;
        _next.parameters[n] += numDiffStep;
        _objective->computeResiduals(_next.parameters, _next.residuals);
        _jacobian.col(n) = (_next.residuals.asEigen() - _current.residuals.asEigen()) / numDiffStep;
        _next.parameters[n] = _current.parameters[n];
    }
    _gradient.deep() = 0.0;
    _hessian.deep() = 0.0;
    if (_objective->hasPrior()) {
        _objective->differentiatePrior(_current.parameters, _gradient, _hessian);
        // objective evaluates P(x); we want -ln P(x) and associated derivatives
        _gradient.asEigen() /= -_current.priorValue;
        _hessian.asEigen() /= -_current.priorValue;
        _hessian.asEigen().selfadjointView<Eigen::Lower>().rankUpdate(_gradient.asEigen(), 1.0);
    }
    if (!_ctrl.noSR1Term) {
        _sr1jtr = _jacobian.adjoint() * _current.residuals.asEigen();
        _gradient.asEigen() += _sr1jtr;
    } else {
        _gradient.asEigen() += _jacobian.adjoint() * _current.residuals.asEigen();
    }
    _hessian.asEigen().selfadjointView<Eigen::Lower>().rankUpdate(_jacobian.adjoint(), 1.0);
}

void Optimizer::removeSR1Term() {
   _hessian.asEigen() -= _sr1b;
}

bool Optimizer::_stepImpl(
    int outerIterCount,
    HistoryRecorder const * recorder,
    afw::table::BaseCatalog * history
) {
    pex::logging::Debug log("meas.modelfit.optimizer.Optimizer");
    _state &= ~int(STATUS);
    if (_gradient.asEigen().lpNorm<Eigen::Infinity>() <= _ctrl.gradientThreshold) {
        log.debug<7>("max(gradient)=%g below threshold; declaring convergence",
                     _gradient.asEigen().lpNorm<Eigen::Infinity>());
        _state |= CONVERGED_GRADZERO;
        return false;
    }
    for (int innerIterCount = 0; innerIterCount < _ctrl.maxInnerIterations; ++innerIterCount) {
        log.debug<10>("Starting inner iteration %d", innerIterCount);
        _state &= ~int(STATUS);
        _next.objectiveValue = 0.0;
        _next.priorValue = 1.0;
        solveTrustRegion(
            _step, _hessian, _gradient, _trustRadius, _ctrl.trustRegionSolverTolerance
        );
        _next.parameters.asEigen() = _current.parameters.asEigen() + _step.asEigen();
        double stepLength = _step.asEigen().norm();
        if (std::isnan(stepLength)) {
            log.debug<7>("NaN encountered in step length");
            _state |= FAILED_NAN;
            return false;
        }
        log.debug<10>("Step has length %g", stepLength);
        if (_objective->hasPrior()) {
            _next.priorValue = _objective->computePrior(_next.parameters);
            _next.objectiveValue = -std::log(_next.priorValue);
            if (_next.priorValue <= 0.0 || std::isnan(_next.objectiveValue)) {
                _next.objectiveValue = std::numeric_limits<Scalar>::infinity();
                log.debug<10>("Rejecting step due to zero prior");
                if (stepLength < _trustRadius) {
                    // Because the step failed due to the prior, we could add an API to the objective
                    // for projecting the step back to a feasible value, and move directly there.
                    // It remains to be seen if that's needed.
                    log.debug<10>("Unconstrained step failed; setting trust radius to step length %g",
                                  stepLength);
                    _trustRadius = stepLength;
                }
                _trustRadius *= _ctrl.trustRegionShrinkFactor;
                log.debug<10>("Decreasing trust radius to %g", _trustRadius);
                _state |= STATUS_STEP_REJECTED | STATUS_TR_DECREASED;
                if (_trustRadius <= _ctrl.minTrustRadiusThreshold) {
                    log.debug<7>("Trust radius %g has dropped below threshold %g; declaring convergence",
                                 _trustRadius, _ctrl.minTrustRadiusThreshold);
                    _state |= CONVERGED_TR_SMALL;
                    return false;
                }
                if (recorder) recorder->apply(outerIterCount, innerIterCount, *history, *this);
                continue;
            }
        }
        _objective->computeResiduals(_next.parameters, _next.residuals);
        _next.objectiveValue += 0.5*_next.residuals.asEigen().squaredNorm();
        double actualChange = _next.objectiveValue - _current.objectiveValue;
        double predictedChange = _step.asEigen().dot(
            _gradient.asEigen() + 0.5*_hessian.asEigen()*_step.asEigen()
        );
        double rho = actualChange / predictedChange;
        if (std::isnan(rho)) {
            log.debug<10>("NaN encountered in rho");
            _state |= FAILED_NAN;
            return false;
        }
        log.debug<10>("Reduction ratio rho=%g; actual=%g, predicted=%g", rho, actualChange, predictedChange);
        if (rho > _ctrl.stepAcceptThreshold && actualChange < 0.0) {
            log.debug<10>("Step accepted; new objective=%g, old was %g",
                          _next.objectiveValue, _current.objectiveValue);
            _state |= STATUS_STEP_ACCEPTED;
            _current.swap(_next);
            if (!_ctrl.noSR1Term) {
                _sr1v = -_sr1jtr;
            }
            _computeDerivatives();
            if (!_ctrl.noSR1Term) {
                _sr1v += _sr1jtr;
                double vs = _sr1v.dot(_step.asEigen());
                if (vs >= (_ctrl.skipSR1UpdateThreshold * _sr1v.norm() * stepLength + 1.0)) {
                    _sr1b.selfadjointView<Eigen::Lower>().rankUpdate(_sr1v, 1.0 / vs);
                }
                _hessian.asEigen() += _sr1b;
            }
            _hessian.asEigen() = _hessian.asEigen().selfadjointView<Eigen::Lower>();
            if (
                rho > _ctrl.trustRegionGrowReductionRatio &&
                (stepLength/_trustRadius) > _ctrl.trustRegionGrowStepFraction
            ) {
                _state |= STATUS_TR_INCREASED;
                _trustRadius *= _ctrl.trustRegionGrowFactor;
                log.debug<10>("Increasing trust radius to %g", _trustRadius);
            } else if (rho < _ctrl.trustRegionShrinkReductionRatio) {
                // even though the step was accepted, our quadratic model
                // of the objective function wasn't very accurate, so we
                // decrease the trust region anyway
                _state |= STATUS_TR_DECREASED;
                _trustRadius *= _ctrl.trustRegionShrinkFactor;
                log.debug<10>("Decreasing trust radius to %g", _trustRadius);
            } else {
                log.debug<10>("Leaving trust radius unchanged at %g", _trustRadius);
                _state |= STATUS_TR_UNCHANGED;
            }
            if (recorder) recorder->apply(outerIterCount, innerIterCount, *history, *this);
            return true;
        }
        _state |= STATUS_STEP_REJECTED;
        log.debug<10>("Step rejected; test objective was %g, current is %g",
                      _next.objectiveValue, _current.objectiveValue);
        if (stepLength < _trustRadius) {
            log.debug<10>("Unconstrained step failed; setting trust radius to step length %g",
                          stepLength);
            _trustRadius = stepLength;
        }
        // we always decrease the trust radius if the step is rejected - otherwise we'll just
        // produce the same step again
        _state |= STATUS_TR_DECREASED;
        _trustRadius *= _ctrl.trustRegionShrinkFactor;
        log.debug<10>("Decreasing trust radius to %g", _trustRadius);
        if (_trustRadius <= _ctrl.minTrustRadiusThreshold) {
            _state |= CONVERGED_TR_SMALL;
            log.debug<7>("Trust radius %g has dropped below threshold %g; declaring convergence",
                         _trustRadius, _ctrl.minTrustRadiusThreshold);
            return false;
        }
        if (recorder) recorder->apply(outerIterCount, innerIterCount, *history, *this);
    }
    log.debug<7>("Max inner iteration number exceeded");
    _state |= FAILED_MAX_INNER_ITERATIONS;
    return false;
}

int Optimizer::_runImpl(HistoryRecorder const * recorder, afw::table::BaseCatalog * history) {
    pex::logging::Debug log("meas.modelfit.optimizer.Optimizer");
    if (recorder) recorder->apply(-1, -1, *history, *this);
    int outerIterCount = 0;
    try {
        for (; outerIterCount < _ctrl.maxOuterIterations; ++outerIterCount) {
            log.debug<10>("Starting outer iteration %d", outerIterCount);
            if (!_stepImpl(outerIterCount, recorder, history)) return outerIterCount;
        }
        _state |= FAILED_MAX_OUTER_ITERATIONS;
        log.debug<7>("Max outer iteration number exceeded");
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
    pex::logging::Debug log("meas.modelfit.optimizer.solveTrustRegion");
    double const r2 = r*r;
    double const r2min = r2 * (1.0 - tolerance) * (1.0 - tolerance);
    double const r2max = r2 * (1.0 + tolerance) * (1.0 + tolerance);
    int const d = g.getSize<0>();
    Eigen::SelfAdjointEigenSolver<Matrix> eigh(F.asEigen());
    double const threshold = ROOT_EPS * eigh.eigenvalues()[d - 1];
    Vector qtg = eigh.eigenvectors().adjoint() * g.asEigen();
    Vector tmp(d);
    double mu = 0.0;
    double xsn = 0.0;
    if (eigh.eigenvalues()[0] >= threshold) {
        log.debug<10>("Starting with full-rank matrix");
        tmp = (eigh.eigenvalues().array().inverse() * qtg.array()).matrix();
        x.asEigen() = -eigh.eigenvectors() * tmp;
        xsn = x.asEigen().squaredNorm();
        if (xsn <= r2max) {
            log.debug<10>("Ending with unconstrained solution");
            // unconstrained solution is within the constraint; no more work to do
            return;
        }
    } else {
        mu = -eigh.eigenvalues()[0] + 2.0*ROOT_EPS*eigh.eigenvalues()[d - 1];
        tmp = ((eigh.eigenvalues().array() + mu).inverse() * qtg.array()).matrix();
        int n = 0;
        while (eigh.eigenvalues()[++n] < threshold);
        log.debug<10>("Starting with %d zero eigenvalue(s) (of %d)", n, d);
        if ((qtg.head(n).array() < ROOT_EPS * g.asEigen().lpNorm<Eigen::Infinity>()).all()) {
            x.asEigen() = -eigh.eigenvectors().rightCols(n) * tmp.tail(n);
            xsn = x.asEigen().squaredNorm();
            if (xsn < r2min) {
                // Nocedal and Wright's "Hard Case", which is actually
                // easier: Q_1^T g is zero (where the columns of Q_1
                // are the eigenvectors that correspond to the
                // smallest eigenvalue \lambda_1), so \mu = -\lambda_1
                // and we can add a multiple of any column of Q_1 to x
                // to get ||x|| == r.  If ||x|| > r, we can find the
                // solution with the usual iteration by increasing \mu.
                double tau = std::sqrt(r*r - x.asEigen().squaredNorm());
                x.asEigen() += tau * eigh.eigenvectors().col(0);
                log.debug<10>("Ending; Q_1^T g == 0, and ||x|| < r");
                return;
            }
            log.debug<10>("Continuing; Q_1^T g == 0, but ||x|| > r");
        } else {
            x.asEigen() = -eigh.eigenvectors() * tmp;
            xsn = x.asEigen().squaredNorm();
            log.debug<10>("Continuing; Q_1^T g != 0, ||x||=%f");
        }
    }
    int nIter = 0;
    while ((xsn < r2min || xsn > r2max) && ++nIter < ITER_MAX) {
        log.debug<10>("Iterating at mu=%f, ||x||=%f, r=%f", mu, std::sqrt(xsn), r);
        mu += xsn*(std::sqrt(xsn) / r - 1.0)
            / (qtg.array().square() / (eigh.eigenvalues().array() + mu).cube()).sum();
        tmp = ((eigh.eigenvalues().array() + mu).inverse() * qtg.array()).matrix();
        x.asEigen() = -eigh.eigenvectors() * tmp;
        xsn = x.asEigen().squaredNorm();
    }
    log.debug<10>("Ending at mu=%f, ||x||=%f, r=%f", mu, std::sqrt(xsn), r);
    return;
}

}}} // namespace lsst::meas::modelfit
