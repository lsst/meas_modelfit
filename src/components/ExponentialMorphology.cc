// -*- lsst-c++ -*-
/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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
 
#include "lsst/meas/multifit/components/ExponentialMorphology.h"
#include "lsst/meas/multifit/components/ExponentialMorphologyProjection.h"
#include "lsst/afw/geom/ellipses/LogShear.h"

namespace components = lsst::meas::multifit::components;

components::ExponentialMorphology::Ptr components::ExponentialMorphology::create(
    Parameter const & flux,
    lsst::afw::geom::ellipses::BaseCore const & ellipse
) { 
    lsst::afw::geom::ellipses::LogShear logShear(ellipse);
    double radius = std::exp(logShear[lsst::afw::geom::ellipses::LogShear::KAPPA]);
    boost::shared_ptr<ParameterVector> linear(new ParameterVector(LINEAR_SIZE));
    boost::shared_ptr<ParameterVector> nonlinear(new ParameterVector(NONLINEAR_SIZE));
    *linear << flux;
    *nonlinear << logShear[GAMMA1], logShear[GAMMA2], radius;
    return ExponentialMorphology::Ptr(new ExponentialMorphology(linear, nonlinear));    
}

lsst::afw::geom::ellipses::Core::Ptr 
components::ExponentialMorphology::computeBoundingEllipseCore() const {  
    ParameterConstIterator params(beginNonlinear());
    return boost::make_shared<lsst::afw::geom::ellipses::LogShear> (
        params[GAMMA1], params[GAMMA2], std::log(params[RADIUS])
    );
}
components::Morphology::Ptr components::ExponentialMorphology::create(
    boost::shared_ptr<ParameterVector const> const & linearParameters,
    boost::shared_ptr<ParameterVector const> const & nonlinearParameters,
    size_t const & start
) const {
    return Morphology::Ptr(
        new ExponentialMorphology(linearParameters, nonlinearParameters, start)
    );
}

components::MorphologyProjection::Ptr components::ExponentialMorphology::makeProjection(
    lsst::afw::geom::Extent2I const & kernelDimensions,
    lsst::afw::geom::AffineTransform const & transform
) const {
    return ExponentialMorphologyProjection::Ptr(new ExponentialMorphologyProjection(
        boost::static_pointer_cast<ExponentialMorphology const>(shared_from_this()),
        kernelDimensions,
        transform
    ));
}

void components::ExponentialMorphology::_handleNonlinearParameterChange() {}
