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
 
#include "lsst/meas/multifit/components/SersicMorphology.h"
#include "lsst/meas/multifit/components/SersicMorphologyProjection.h"
#include "lsst/afw/geom/ellipses/LogShear.h"
#include "lsst/pex/exceptions/Runtime.h"

namespace components = lsst::meas::multifit::components;

lsst::meas::multifit::Cache::ConstPtr components::SersicMorphology::_cache;

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
    lsst::afw::geom::AffineTransform::ConstPtr const & transform
) const {
    return SersicMorphologyProjection::Ptr(new SersicMorphologyProjection(
        boost::static_pointer_cast<SersicMorphology const>(shared_from_this()),
        kernelDimensions,
        transform
    ));
}

void components::SersicMorphology::_handleNonlinearParameterChange() {
    if(!_cache) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Must call SersicMorphology::setSersicCache before creating and using SersicMorphology"
        );
    }

    ParameterMap errors;
    int offset = getNonlinearParameterOffset();

    lsst::afw::geom::ellipses::Distortion dist(*computeBoundingEllipseCore());
    if(dist.getE() > 0.95) {
        dist.setE(0.95);
        lsst::afw::geom::ellipses::LogShear ls(dist);
        Parameter gamma1 = ls[lsst::afw::geom::ellipses::LogShear::GAMMA1];
        Parameter gamma2 = ls[lsst::afw::geom::ellipses::LogShear::GAMMA2];
        errors[offset + GAMMA1] = gamma1;
        errors[offset + GAMMA2] = gamma2;
    }

    Parameter sersic = getSersicIndex();
    lsst::afw::geom::BoxD bounds = _cache->getParameterBounds();
    try {
        _indexFunctor = _cache->getRowFunctor(sersic);
    } 
    catch(lsst::pex::exceptions::InvalidParameterException &e) {
        sersic = std::min(std::max(sersic, bounds.getMinY()), bounds.getMaxY()); 
        errors[offset + SERSIC_INDEX] = sersic;
    }

    if(errors.size() > 0) {
        throw LSST_EXCEPT(
            lsst::meas::multifit::ParameterRangeException,
            "SersicMorphology nonlinear parameters out of range", errors
        );
    }
}
