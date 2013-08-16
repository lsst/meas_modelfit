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

#include "Eigen/LU"

#include "lsst/pex/exceptions.h"
#include "lsst/meas/multifit/priors.h"
#include "lsst/meas/multifit/integrals.h"

namespace lsst { namespace meas { namespace multifit {

samples::Scalar FlatPrior::apply(LogGaussian const & likelihood, samples::Vector const & parameters) const {
    return integrateGaussian(likelihood.grad, likelihood.fisher)
        + std::log(_maxRadius * _maxEllipticity * _maxEllipticity * 2 * M_PI);
}

FlatPrior::FlatPrior(double maxRadius, double maxEllipticity) :
    Prior(ParameterDefinition::lookup("SeparableReducedShearTraceRadius")),
    _maxRadius(maxRadius),
    _maxEllipticity(maxEllipticity)
{}

}}} // namespace lsst::meas::multifit
