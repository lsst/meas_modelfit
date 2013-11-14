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

#include <Eigen/Eigenvalues>

#include "ndarray/eigen.h"

#define LSST_MAX_DEBUG 10
#include "lsst/pex/logging/Debug.h"
#include "lsst/utils/ieee.h"

#include "lsst/meas/multifit/optimizer.h"
#include "lsst/meas/multifit/Likelihood.h"
#include "lsst/meas/multifit/priors.h"

namespace lsst { namespace meas { namespace multifit {

void solveTrustRegion(
    ndarray::Array<Scalar,1,1> const & x,
    ndarray::Array<Scalar const,2,1> const & F,
    ndarray::Array<Scalar const,1,1> const & g,
    double r, double tolerance
) {
    static double const ROOT_EPS = std::sqrt(std::numeric_limits<double>::epsilon());
    static int const ITER_MAX = 10;
    pex::logging::Debug log("meas.multifit.optimizer.solveTrustRegion");
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
        log.debug<7>("Starting with full-rank matrix");
        tmp = (eigh.eigenvalues().array().inverse() * qtg.array()).matrix();
        x.asEigen() = -eigh.eigenvectors() * tmp;
        xsn = x.asEigen().squaredNorm();
        if (xsn <= r2max) {
            log.debug<7>("Ending with unconstrained solution");
            // unconstrained solution is within the constraint; no more work to do
            return;
        }
    } else {
        mu = -eigh.eigenvalues()[0] + 2.0*ROOT_EPS*eigh.eigenvalues()[d - 1];
        tmp = ((eigh.eigenvalues().array() + mu).inverse() * qtg.array()).matrix();
        int n = 0;
        while (eigh.eigenvalues()[++n] < threshold);
        log.debug<7>("Starting with %d zero eigenvalue(s) (of %d)", n, d);
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
                log.debug<7>("Ending; Q_1^T g == 0, and ||x|| < r");
                return;
            }
            log.debug<7>("Continuing; Q_1^T g == 0, but ||x|| > r");
        } else {
            x.asEigen() = -eigh.eigenvectors() * tmp;
            xsn = x.asEigen().squaredNorm();
            log.debug<7>("Continuing; Q_1^T g != 0");
        }
    }
    int nIter = 0;
    while ((xsn < r2min || xsn > r2max) && ++nIter < ITER_MAX) {
        log.debug<7>("Iterating at mu=%f, ||x||=%f, r=%f", mu, std::sqrt(xsn), r);
        mu += xsn*(std::sqrt(xsn) / r - 1.0)
            / (qtg.array().square() / (eigh.eigenvalues().array() + mu).cube()).sum();
        tmp = ((eigh.eigenvalues().array() + mu).inverse() * qtg.array()).matrix();
        x.asEigen() = -eigh.eigenvectors() * tmp;
        xsn = x.asEigen().squaredNorm();
    }
    log.debug<7>("Ending at mu=%f, ||x||=%f, r=%f", mu, std::sqrt(xsn), r);
    return;
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

    virtual bool hasPrior() const { return _prior; }

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
    return boost::make_shared<LikelihoodOptimizerObjective>(likelihood, prior);
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
    pex::logging::Debug log("meas.multifit.optimizer.Optimizer");
    if (parameters.getSize<0>() != _objective->parameterSize) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthErrorException,
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
    log.debug<6>("Initial objective value is %g", _current.objectiveValue);
    _sr1b.setZero();
    _computeDerivatives();
    _hessian.asEigen() = _hessian.asEigen().selfadjointView<Eigen::Lower>();
}

void Optimizer::_computeDerivatives() {
    _jacobian.setZero();
    _next.parameters.deep() = _current.parameters;
    for (int n = 0; n < _objective->parameterSize; ++n) {
        double numDiffStep = _ctrl.numDiffRelStep * _next.parameters[n] + _ctrl.numDiffAbsStep;
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

bool Optimizer::step() {
    pex::logging::Debug log("meas.multifit.optimizer.Optimizer");
    _state &= ~int(STATUS);
    if (_gradient.asEigen().lpNorm<Eigen::Infinity>() <= _ctrl.gradientThreshold) {
        log.debug<6>("max(gradient)=%g below threshold; declaring convergence",
                     _gradient.asEigen().lpNorm<Eigen::Infinity>());
        _state |= CONVERGED_GRADZERO;
        return false;
    }
    for (int iterCount = 0; iterCount < _ctrl.maxInnerIterations; ++iterCount) {
        log.debug<6>("Starting inner iteration %d", iterCount);
        _state &= ~int(STATUS);
        _next.objectiveValue = 0.0;
        _next.priorValue = 1.0;
        solveTrustRegion(
            _step, _hessian, _gradient, _trustRadius, _ctrl.trustRegionSolverTolerance
        );
        _next.parameters.asEigen() = _current.parameters.asEigen() + _step.asEigen();
        double stepLength = _step.asEigen().norm();
        if (utils::isnan(stepLength)) {
            log.debug<6>("NaN encountered in step length");
            std::cerr << "gradient:\n" << _gradient << std::endl;
            std::cerr << "hessian:\n" << _hessian << std::endl;
            _state |= FAILED_NAN;
            return false;
        }
        log.debug<6>("Step has length %g", stepLength);
        if (_objective->hasPrior()) {
            _next.priorValue = _objective->computePrior(_next.parameters);
            if (_next.priorValue <= 0.0) {
                log.debug<6>("Rejecting step due to zero prior");
                _trustRadius *= _ctrl.trustRegionShrinkFactor;
                log.debug<6>("Decreasing trust radius to %g", _trustRadius);
                _state |= STATUS_STEP_REJECTED | STATUS_TR_DECREASED;
                if (_trustRadius <= _ctrl.minTrustRadiusThreshold) {
                    log.debug<6>("Trust radius %g has dropped below threshold %g; declaring convergence",
                                 _trustRadius, _ctrl.minTrustRadiusThreshold);
                    _state |= CONVERGED_TR_SMALL;
                    return false;
                }
                continue;
            }
            _next.objectiveValue = -std::log(_next.priorValue);
        }
        _objective->computeResiduals(_next.parameters, _next.residuals);
        _next.objectiveValue += 0.5*_next.residuals.asEigen().squaredNorm();
        double actualChange = _next.objectiveValue - _current.objectiveValue;
        double predictedChange = _step.asEigen().dot(
            _gradient.asEigen() + 0.5*_hessian.asEigen()*_step.asEigen()
        );
        double rho = actualChange / predictedChange;
        if (utils::isnan(rho)) {
            log.debug<6>("NaN encountered in rho");
            _state |= FAILED_NAN;
            return false;
        }
        log.debug<6>("Reduction ratio rho=%g; actual=%g, predicted=%g", rho, actualChange, predictedChange);
        if (rho > _ctrl.stepAcceptThreshold) {
            log.debug<6>("Step accepted");
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
                log.debug<6>("Increasing trust radius to %g", _trustRadius);
            } else {
                log.debug<6>("Leaving trust radius unchanged at %g", _trustRadius);
                _state |= STATUS_TR_UNCHANGED;
            }
            return true;
        }
        log.debug<6>("Step rejected");
        if (rho < _ctrl.trustRegionShrinkReductionRatio) {
            _state |= STATUS_TR_DECREASED;
            _trustRadius *= _ctrl.trustRegionShrinkFactor;
            log.debug<6>("Decreasing trust radius to %g", _trustRadius);
            if (_trustRadius <= _ctrl.minTrustRadiusThreshold) {
                _state |= CONVERGED_TR_SMALL;
                log.debug<6>("Trust radius %g has dropped below threshold %g; declaring convergence",
                             _trustRadius, _ctrl.minTrustRadiusThreshold);
                return false;
            }
        } else {
            log.debug<6>("Leaving trust radius unchanged at %g", _trustRadius);
            _state |= STATUS_TR_UNCHANGED;
        }
    }
    log.debug<6>("Max inner iteration number exceeded");
    _state |= FAILED_MAX_INNER_ITERATIONS;
    return false;
}

int Optimizer::run() {
    pex::logging::Debug log("meas.multifit.optimizer.Optimizer");
    int iterCount = 0;
    try {
        for (; iterCount < _ctrl.maxOuterIterations; ++iterCount) {
            log.debug<6>("Starting outer iteration %d", iterCount);
            if (!step()) return iterCount;
        }
        _state |= FAILED_MAX_OUTER_ITERATIONS;
    } catch (...) {
        _state |= FAILED_EXCEPTION;
    }
    return iterCount;
}



}}} // namespace lsst::meas::multifit
