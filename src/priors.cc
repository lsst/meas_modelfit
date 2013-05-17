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

#include "boost/math/special_functions/erf.hpp"

#include "lsst/pex/exceptions.h"
#include "lsst/meas/multifit/priors.h"

namespace lsst { namespace meas { namespace multifit {

namespace {

// a 1-d integrals we use repeatedly below
// @f$\int_0^{\infty} \alpha^n e^{-(\alpha-\mu)^2 f / 2} d\alpha
double integral(int n, double mu, double f) {
    double i0 = std::sqrt(0.5 * M_PI / f) * (1.0 + boost::math::erf(std::sqrt(0.5 * f) * mu));
    if (n == 0) {
        return i0;
    } else {
        double t1 = std::exp(-0.5 * f * mu * mu) / f;
        if (n == 1) {
            return t1 + i0 * mu;
        }
        if (n == 2) {
            return t1 * mu + i0 * (mu * mu + 1.0 / f);
        }
    }
    throw LSST_EXCEPT(
        pex::exceptions::LogicErrorException,
        "Moment must be <= 2"
    );
}

} // anonymous

double SingleComponentPrior::apply(LogGaussian const & likelihood, Vector const & parameters) const {
    if (likelihood.mu.size() != 2) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "SingleComponentPrior is only valid for two-component models"
        );
    }
    return std::exp(-0.5 * likelihood.r) * (
        _beta * integral(0, likelihood.mu[0], likelihood.fisher(0,0))
        + (1.0 - _beta) * integral(0, likelihood.mu[1], likelihood.fisher(1,1))
    );
}

double SingleComponentPrior::computeFluxExpectation(
    LogGaussian const & likelihood, Vector const & parameters
) const {
    if (likelihood.mu.size() != 2) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "SingleComponentPrior is only valid for two-component models"
        );
    }
    return std::exp(-0.5 * likelihood.r) * (
        _beta * integral(1, likelihood.mu[0], likelihood.fisher(0,0))
        + (1.0 - _beta) * integral(1, likelihood.mu[1], likelihood.fisher(1,1))
    );
}

double SingleComponentPrior::computeSquaredFluxExpectation(
    LogGaussian const & likelihood, Vector const & parameters
) const {
    if (likelihood.mu.size() != 2) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "SingleComponentPrior is only valid for two-component models"
        );
    }
    return std::exp(-0.5 * likelihood.r) * (
        _beta * integral(2, likelihood.mu[0], likelihood.fisher(0,0))
        + (1.0 - _beta) * integral(2, likelihood.mu[1], likelihood.fisher(1,1))
    );
}

Vector SingleComponentPrior::computeFractionExpectation(
    LogGaussian const & likelihood, Vector const & parameters
) const {
    if (likelihood.mu.size() != 2) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "SingleComponentPrior is only valid for two-component models"
        );
    }
    Vector result(2);
    result[0] = _beta * integral(0, likelihood.mu[0], likelihood.fisher(0,0));
    result[1] = (1 - _beta) * integral(0, likelihood.mu[1], likelihood.fisher(1,1));
    return result;
}

}}} // namespace lsst::meas::multifit
