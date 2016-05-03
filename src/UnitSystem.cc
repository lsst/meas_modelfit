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

#include "boost/format.hpp"
#include "boost/make_shared.hpp"

#include "lsst/meas/modelfit/UnitSystem.h"

namespace lsst { namespace meas { namespace modelfit {

UnitSystem::UnitSystem(afw::coord::Coord const & position, Scalar mag) {
    double cdelt = (1.0*lsst::afw::geom::arcseconds).asDegrees();
    wcs = afw::image::makeWcs(position, afw::geom::Point2D(0.0, 0.0), cdelt, 0.0, 0.0, cdelt);
    PTR(afw::image::Calib) calib_ = std::make_shared<afw::image::Calib>();
    calib_->setFluxMag0(std::pow(10.0, mag/2.5));
    calib = calib_;
}

LocalUnitTransform::LocalUnitTransform(
    afw::coord::Coord const & position,
    UnitSystem const & source,
    UnitSystem const & destination
) :
    geometric(
        afw::image::XYTransformFromWcsPair(
            destination.wcs,
            source.wcs
        ).linearizeForwardTransform(source.wcs->skyToPixel(position))
    ),
    flux(destination.calib->getFluxMag0().first / source.calib->getFluxMag0().first),
    sb(flux / geometric.getLinear().computeDeterminant())
{}

}}} // namespace lsst::meas::modelfit
