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

#include "lsst/meas/multifit/optimizer.h"

namespace lsst { namespace meas { namespace multifit {

void solveTrustRegion(
    ndarray::Array<double,1,1> const & x,
    ndarray::Array<double const,2,1> const & F,
    ndarray::Array<double const,1,1> const & g,
    double r, double tolerance
) {
    static double const ROOT_EPS = std::sqrt(std::numeric_limits<double>::epsilon());
    static int const ITER_MAX = 10;
    pex::logging::Debug log("meas.multifit.optimizer.solveTrustRegion");
    double const r2 = r*r;
    double const r2min = r2 * (1.0 - tolerance) * (1.0 - tolerance);
    double const r2max = r2 * (1.0 + tolerance) * (1.0 + tolerance);
    int const d = g.getSize<0>();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigh(F.asEigen());
    double const threshold = ROOT_EPS * eigh.eigenvalues()[d - 1];
    Eigen::VectorXd qtg = eigh.eigenvectors().adjoint() * g.asEigen();
    Eigen::VectorXd tmp(d);
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
        log.debug<7>("Starting with %d zero eigenvalue(s)", n);
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


double PosteriorOptimizerObjective::computePrior(ndarray::Array<double const,1,1> const & parameters) const {
    return 1.0;
}

void PosteriorOptimizerObjective::differentiatePrior(
    ndarray::Array<double const,1,1> const & parameters,
    ndarray::Array<double,1,1> const & gradient,
    ndarray::Array<double,2,1> const & hessian
) const {
    gradient.deep() = 0.0;
    hessian.deep() = 0.0;
}

PosteriorOptimizerIterationData::PosteriorOptimizerIterationData(int dataSize, int parameterSize) :
    objectiveValue(0.0), priorValue(0.0),
    parameters(ndarray::allocate(parameterSize)),
    residuals(ndarray::allocate(dataSize))
{}

void PosteriorOptimizerIterationData::swap(PosteriorOptimizerIterationData & other) {
    std::swap(objectiveValue, other.objectiveValue);
    std::swap(priorValue, other.priorValue);
    parameters.swap(other.parameters);
    residuals.swap(other.residuals);
}


PosteriorOptimizer::PosteriorOptimizer(
    PTR(Objective const) objective,
    ndarray::Array<double const,1,1> const & parameters,
    Control const & ctrl
) :
    _state(0x0),
    _objective(objective),
    _ctrl(ctrl),
    _trustRadius(std::numeric_limits<double>::infinity()),
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
    _sr1b.setZero();
    _computeDerivatives();
    _hessian.asEigen() = _hessian.asEigen().selfadjointView<Eigen::Lower>();
}

void PosteriorOptimizer::_computeDerivatives() {
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

bool PosteriorOptimizer::step() {
    for (int iterCount = 0; iterCount < _ctrl.maxInnerIterations; ++iterCount) {
        _state &= ~int(STATUS);
        _state |= FAILED_EXCEPTION; // set in advance in case we throw
        solveTrustRegion(
            _step, _hessian, _gradient, _trustRadius, _ctrl.trustRegionSolverTolerance
        );
        _next.parameters.asEigen() = _current.parameters.asEigen() + _step.asEigen();
        double stepLength = _step.asEigen().norm();
        if (_trustRadius == std::numeric_limits<double>::infinity()) {
            _trustRadius = stepLength;
        }
        if (_objective->hasPrior()) {
            _next.priorValue = _objective->computePrior(_next.parameters);
            if (_next.priorValue <= 0.0) {
                _trustRadius *= _ctrl.trustRegionShrinkFactor;
                _state |= STATUS_STEP_REJECTED | STATUS_TR_DECREASED;
                if (_trustRadius <= _ctrl.minTrustRadiusThreshold) {
                    _state |= CONVERGED_TR_SMALL;
                    return true;
                }
                continue;
            }
            _next.objectiveValue = -std::log(_next.priorValue);
        }
        _objective->computeResiduals(_next.parameters, _next.residuals);
        _next.objectiveValue += 0.5*_next.residuals.asEigen().squaredNorm();
        double actualReduction = _current.objectiveValue - _next.objectiveValue;
        double predictedReduction = -_step.asEigen().dot(
            _gradient.asEigen() + 0.5*_hessian.asEigen()*_step.asEigen()
        );
        double rho = actualReduction / predictedReduction;
        if (rho > _ctrl.stepAcceptThreshold) {
            _state |= STATUS_STEP_ACCEPTED;
            _current.swap(_next);
            if (!_ctrl.noSR1Term) {
                _sr1v = -_sr1jtr;
            }
            _computeDerivatives();
            if (!_ctrl.noSR1Term) {
                _sr1v += _sr1jtr;
                double vs = _sr1v.dot(_step.asEigen());
                if (vs >= _ctrl.skipSR1UpdateThreshold * _sr1v.norm() * stepLength) {
                    _sr1b.selfadjointView<Eigen::Lower>().rankUpdate(_sr1v, 1.0 / vs);
                }
                _hessian.asEigen() += _sr1b;
            }
            _hessian.asEigen() = _hessian.asEigen().selfadjointView<Eigen::Lower>();
            if (_gradient.asEigen().lpNorm<Eigen::Infinity>() <= _ctrl.gradientThreshold) {
                _state |= CONVERGED_GRADZERO;
                return true;
            }
            return false;
        }
        if (rho > _ctrl.trustRegionGrowReductionRatio && rho > _ctrl.trustRegionGrowStepFraction) {
            _state |= STATUS_TR_INCREASED;
            _trustRadius *= _ctrl.trustRegionGrowFactor;
        } else if (
            rho > _ctrl.trustRegionShrinkMinReductionRatio
            && rho < _ctrl.trustRegionShrinkMaxReductionRatio
        ) {
            _state |= STATUS_TR_DECREASED;
            _trustRadius *= _ctrl.trustRegionShrinkFactor;
            if (_trustRadius <= _ctrl.minTrustRadiusThreshold) {
                _state |= CONVERGED_TR_SMALL;
                return true;
            }
        } else {
            _state |= STATUS_TR_UNCHANGED;
        }
        _state &= ~int(FAILED_EXCEPTION); // clear the exception flag if we got here
    }
    _state |= FAILED_MAX_INNER_ITERATIONS;
    return true;
}

}}} // namespace lsst::meas::multifit
