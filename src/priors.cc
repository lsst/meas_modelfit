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

#include "lsst/meas/multifit/priors.h"

namespace lsst { namespace meas { namespace multifit {

double SingleComponentPrior::apply(LogGaussian const & likelihood, Vector const & parameters) const {
    if (likelihood.mu.size() != 2) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "SingleComponentPrior is only valid for two-component models"
        );
    }
    double m1 = _beta * std::sqrt(0.5 * M_PI / likelihood.fisher(0,0))
        * (1.0 + boost::math::erf(std::sqrt(0.5 * likelihood.fisher(0,0)) * likelihood.mu[0]));
    double m2 = (1.0 - _beta) * std::sqrt(0.5 * M_PI / likelihood.fisher(1,1))
        * (1.0 + boost::math::erf(std::sqrt(0.5 * likelihood.fisher(1,1)) * likelihood.mu[1]));
    return std::exp(-0.5 * likelihood.r) * (m1 + m2);
}

}}} // namespace lsst::meas::multifit
