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

#include <boost/format.hpp>
#include "lsst/meas/multifit/components/SersicMorphology.h"
#include "lsst/meas/multifit/components/SersicMorphologyProjection.h"
#include "lsst/afw/geom/ellipses/LogShear.h"
#include "lsst/pex/exceptions/Runtime.h"

namespace components = lsst::meas::multifit::components;

lsst::meas::multifit::SersicCache::ConstPtr components::SersicMorphology::_cache;

lsst::afw::geom::ellipses::Core::Ptr 
components::SersicMorphology::computeBoundingEllipseCore() const {  
    ParameterConstIterator params(beginNonlinear());
    return boost::make_shared<lsst::afw::geom::ellipses::LogShear> (
        params[GAMMA1], params[GAMMA2], params[KAPPA]
    );
}
components::Morphology::Ptr components::SersicMorphology::create(
    boost::shared_ptr<ParameterVector const> const & linearParameters,
    boost::shared_ptr<ParameterVector const> const & nonlinearParameters,
    size_t const & start
) const {
    return Morphology::Ptr(
        new SersicMorphology(linearParameters, nonlinearParameters, start)
    );
}

components::MorphologyProjection::Ptr components::SersicMorphology::makeProjection(
    lsst::afw::geom::Extent2I const & kernelDimensions,
    lsst::afw::geom::AffineTransform const & transform
) const {
    return SersicMorphologyProjection::Ptr(new SersicMorphologyProjection(
        boost::static_pointer_cast<SersicMorphology const>(shared_from_this()),
        kernelDimensions,
        transform
    ));
}

void components::SersicMorphology::checkCache() {
    if(!_cache) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Must call SersicMorphology::setSersicCache before creating and using SersicMorphology"
        );
    }
}

void components::SersicMorphology::_handleNonlinearParameterChange() {
    checkCache();
    Parameter n = getSersicIndex();
/*
    ParameterMap errors;
    std::string errStr="";
    int offset = getNonlinearParameterOffset();
    lsst::afw::geom::ellipses::Distortion dist(*computeBoundingEllipseCore());
    if(dist.getE() > 0.95) {
        errStr += (boost::format(" Distortion %.3f is too large.")%dist.getE()).str();
        dist.setE(0.95);
        lsst::afw::geom::ellipses::LogShear ls(dist);
        Parameter gamma1 = ls[lsst::afw::geom::ellipses::LogShear::GAMMA1];
        Parameter gamma2 = ls[lsst::afw::geom::ellipses::LogShear::GAMMA2];
        errors[offset + GAMMA1] = gamma1;
        errors[offset + GAMMA2] = gamma2;
    }

    if(n < _cache->getSersicMin() || n >= _cache->getSersicMax()) {
        errStr += (boost::format(" Sersic Index %.3f out of bounds [%.3f, %.3f)")%
                n % _cache->getSersicMin() % _cache->getSersicMax()).str();
        n = std::min(std::max(n, _cache->getSersicMin()), _cache->getSersicMax()); 
        errors[offset+SERSIC_INDEX] = n;
    }

    if(errors.size() > 0) {
        throw LSST_EXCEPT(
            lsst::meas::multifit::ParameterRangeException,
            "SersicMorphology nonlinear parameters out of range."+errStr, errors
        );
    }
*/
    _interpolator = _cache->getInterpolator(n);
    _derivativeInterpolator = _cache->getDerivativeInterpolator(n);
}
