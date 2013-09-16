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

#define LSST_MAX_DEBUG 10
#include "lsst/pex/logging/Debug.h"

#include "lsst/meas/multifit/optimizer.h"

namespace lsst { namespace meas { namespace multifit {

Eigen::VectorXd solveTrustRegion(
    Eigen::MatrixXd const & F, Eigen::VectorXd const & g, double r, double tolerance
) {
    static double const ROOT_EPS = std::sqrt(std::numeric_limits<double>::epsilon());
    static int const ITER_MAX = 10;
    pex::logging::Debug log("meas.multifit.optimizer.solveTrustRegion");
    double const r2 = r*r;
    double const r2min = r2 * (1.0 - tolerance) * (1.0 - tolerance);
    double const r2max = r2 * (1.0 + tolerance) * (1.0 + tolerance);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigh(F);
    double const threshold = ROOT_EPS * eigh.eigenvalues()[g.size() - 1];
    Eigen::VectorXd qtg = eigh.eigenvectors().adjoint() * g;
    Eigen::VectorXd tmp(g.size());
    Eigen::VectorXd x(g.size());
    double mu = 0.0;
    double xsn = 0.0;
    if (eigh.eigenvalues()[0] >= threshold) {
        log.debug<7>("Starting with full-rank matrix");
        tmp = (eigh.eigenvalues().array().inverse() * qtg.array()).matrix();
        x = -eigh.eigenvectors() * tmp;
        xsn = x.squaredNorm();
        if (xsn <= r2max) {
            log.debug<7>("Ending with unconstrained solution");
            // unconstrained solution is within the constraint; no more work to do
            return x;
        }
    } else {
        mu = -eigh.eigenvalues()[0] + 2.0*ROOT_EPS*eigh.eigenvalues()[g.size() - 1];
        tmp = ((eigh.eigenvalues().array() + mu).inverse() * qtg.array()).matrix();
        int n = 0;
        while (eigh.eigenvalues()[++n] < threshold);
        log.debug<7>("Starting with %d zero eigenvalue(s)", n);
        if ((qtg.head(n).array() < ROOT_EPS * g.lpNorm<Eigen::Infinity>()).all()) {
            x = -eigh.eigenvectors().rightCols(n) * tmp.tail(n);
            xsn = x.squaredNorm();
            if (xsn < r2min) {
                // Nocedal and Wright's "Hard Case", which is actually
                // easier: Q_1^T g is zero (where the columns of Q_1
                // are the eigenvectors that correspond to the
                // smallest eigenvalue \lambda_1), so \mu = -\lambda_1
                // and we can add a multiple of any column of Q_1 to x
                // to get ||x|| == r.  If ||x|| > r, we can find the
                // solution with the usual iteration by increasing \mu.
                double tau = std::sqrt(r*r - x.squaredNorm());
                x += tau * eigh.eigenvectors().col(0);
                log.debug<7>("Ending; Q_1^T g == 0, and ||x|| < r");
                return x;
            }
            log.debug<7>("Continuing; Q_1^T g == 0, but ||x|| > r");
        } else {
            x = -eigh.eigenvectors() * tmp;
            xsn = x.squaredNorm();
            log.debug<7>("Continuing; Q_1^T g != 0");
        }
    }
    int nIter = 0;
    while ((xsn < r2min || xsn > r2max) && ++nIter < ITER_MAX) {
        log.debug<7>("Iterating at mu=%f, ||x||=%f, r=%f", mu, std::sqrt(xsn), r);
        mu += xsn*(std::sqrt(xsn) / r - 1.0)
            / (qtg.array().square() / (eigh.eigenvalues().array() + mu).cube()).sum();
        tmp = ((eigh.eigenvalues().array() + mu).inverse() * qtg.array()).matrix();
        x = -eigh.eigenvectors() * tmp;
        xsn = x.squaredNorm();
    }
    log.debug<7>("Ending at mu=%f, ||x||=%f, r=%f", mu, std::sqrt(xsn), r);
    return x;
}


PosteriorOptimizer::PosteriorOptimizer(
    PTR(Objective const) objective,
    Eigen::VectorXd const & parameters,
    Control const & ctrl
) :
    _objective(objective),
    _ctrl(ctrl),
    _trustRadius(std::numeric_limits<double>::infinity()),
    _current(objective->dataSize, objective->parameterSize),
    _next(objective->dataSize, objective->parameterSize),
    _step(objective->parameterSize),
    _jacobian(objective->dataSize, objective->parameterSize),
    _gradient(objective->parameterSize),
    _hessian(objective->parameterSize, objective->parameterSize),
    _sr1b(objective->parameterSize, objective->parameterSize)
{
    if (parameters.size() != _objective->parameterSize) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthErrorException,
            (boost::format("Parameter vector size (%d) does not match objective (%d)")
             % parameters.size() % _objective->parameterSize).str()
        );
    }
    _current.parameters = parameters;
    _next.parameters = parameters;
    _objective->computeResiduals(_current.parameters, _current.residuals);
    _current.objectiveValue = 0.5*_current.residuals.squaredNorm();
    if (_objective->hasPrior()) {
        _current.priorValue = _objective->computePrior(_current.parameters);
        _current.objectiveValue -= std::log(_current.priorValue);
    }
    _sr1b.setZero();
    _computeDerivatives();
    _hessian = _hessian.selfadjointView<Eigen::Lower>();
}

void PosteriorOptimizer::_computeDerivatives() {
    _jacobian.setZero();
    for (int n = 0; n < _objective->parameterSize; ++n) {
        double numDiffStep = _ctrl.numDiffRelStep * _next.parameters[n] + _ctrl.numDiffAbsStep;
        _next.parameters = _current.parameters;
        _next.parameters[n] += numDiffStep;
        _objective->computeResiduals(_next.parameters, _next.residuals);
        _jacobian.col(n) = (_next.residuals - _current.residuals) / numDiffStep;
    }
    _gradient.setZero();
    _hessian.setZero();
    if (_objective->hasPrior()) {
        _objective->differentiatePrior(_current.parameters, _gradient, _hessian);
        // objective evaluates P(x); we want -ln P(x) and associated derivatives
        _gradient /= -_current.priorValue;
        _hessian /= -_current.priorValue;
        _hessian.selfadjointView<Eigen::Lower>().rankUpdate(_gradient, 1.0);
    }
    if (!_ctrl.noSR1Term) {
        _sr1jtr = _jacobian.adjoint() * _current.residuals;
        _gradient += _sr1jtr;
    } else {
        _gradient += _jacobian.adjoint() * _current.residuals;
    }
    _hessian.selfadjointView<Eigen::Lower>().rankUpdate(_jacobian.adjoint(), 1.0);
}

bool PosteriorOptimizer::step() {
    for (int iterCount = 0; iterCount < _ctrl.maxInnerIterations; ++iterCount) {
        _next.parameters = solveTrustRegion(
            _hessian, _gradient, _trustRadius, _ctrl.trustRegionSolverTolerance
        );
        _step = _next.parameters - _current.parameters;
        double stepLength = _step.norm();
        if (_trustRadius == std::numeric_limits<double>::infinity()) {
            _trustRadius = stepLength;
        }
        if (_objective->hasPrior()) {
            _next.priorValue = _objective->computePrior(_next.parameters);
            if (_next.priorValue <= 0.0) {
                _trustRadius *= 0.5;
                continue;
            }
            _next.objectiveValue = -std::log(_next.priorValue);
        }
        _objective->computeResiduals(_next.parameters, _next.residuals);
        _next.objectiveValue += 0.5*_next.residuals.squaredNorm();
        double actualReduction = _current.objectiveValue - _next.objectiveValue;
        double predictedReduction = -(_gradient + 0.5*_hessian*_step).dot(_step);
        double rho = actualReduction / predictedReduction;
        if (rho > _ctrl.stepAcceptThreshold) {
            _current.swap(_next);
            if (!_ctrl.noSR1Term) {
                _sr1v = -_sr1jtr;
            }
            _computeDerivatives();
            if (!_ctrl.noSR1Term) {
                _sr1v += _sr1jtr;
                double vs = _sr1v.dot(_step);
                if (vs >= _ctrl.skipSR1UpdateThreshold * _sr1v.norm() * stepLength) {
                    _sr1b.selfadjointView<Eigen::Lower>().rankUpdate(_sr1v, 1.0 / vs);
                }
                _hessian += _sr1b;
            }
            _hessian = _hessian.selfadjointView<Eigen::Lower>();
            // check convergence
            return true;
        }
        if (rho > _ctrl.trustRegionGrowReductionRatio && rho > _ctrl.trustRegionGrowStepFraction) {
            _trustRadius *= _ctrl.trustRegionGrowFactor;
        }
        if (
            rho > _ctrl.trustRegionShrinkMinReductionRatio && rho < _ctrl.trustRegionShrinkMaxReductionRatio
        ) {
            _trustRadius *= _ctrl.trustRegionShrinkFactor;
        }
    }
    return true;
}

}}} // namespace lsst::meas::multifit
