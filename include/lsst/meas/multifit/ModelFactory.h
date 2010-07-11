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
 
#ifndef LSST_MEAS_MULTIFIT_MODEL_FACTORY_H
#define LSST_MEAS_MULTIFIT_MODEL_FACTORY_H

#include "lsst/afw/geom/Point.h"
#include "lsst/afw/coord/Coord.h"
#include "lsst/meas/multifit/ComponentModel.h"
#include "lsst/meas/multifit/components/Astrometry.h"
#include "lsst/meas/multifit/components/PointSourceMorphology.h"
#include "lsst/meas/multifit/components/SersicMorphology.h"

namespace lsst{
namespace meas{
namespace multifit {

class ModelFactory {
public:
    static Model::Ptr createPointSourceModel(
        Parameter const & flux,
        lsst::afw::geom::Point2D const & centroid
    ) {
        components::Morphology::Ptr morphology = 
            components::PointSourceMorphology::create(flux);
        components::Astrometry::Ptr astrometry(
            new components::Astrometry(centroid)
        );
        return ComponentModel::create(astrometry, morphology);
    }
    static Model::Ptr createPointSourceModel(
        Parameter const & flux,
        lsst::afw::coord::Coord const & coord
    ) {
        return createPointSourceModel(
            flux,
            coord.getPosition(lsst::afw::coord::DEGREES)
        );
    }

    static Model::Ptr createSersicModel(
        Parameter const & flux,
        lsst::afw::geom::Point2D const & centroid,
        lsst::afw::geom::ellipses::Core const & ellipse,
        Parameter const & sersicIndex
    ) {
        components::Morphology::Ptr morphology = 
            components::SersicMorphology::create(flux, ellipse, sersicIndex);
        components::Astrometry::Ptr astrometry(
            new components::Astrometry(centroid)
        );
        return ComponentModel::create(astrometry, morphology); 
    }
    
    static Model::Ptr createSersicModel(
        Parameter const & flux,
        lsst::afw::coord::Coord const & coord,
        lsst::afw::geom::ellipses::Core const & ellipse,
        Parameter const & sersicIndex
    ) {
        return createSersicModel(
            flux, 
            coord.getPosition(lsst::afw::coord::DEGREES),
            ellipse, 
            sersicIndex
        );
    }
};

}}}

#endif
