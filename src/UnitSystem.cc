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

#include <cmath>

#include "boost/format.hpp"
#include <memory>

#include "lsst/afw/geom/transformFactory.h"
#include "lsst/meas/modelfit/UnitSystem.h"

namespace lsst {
namespace meas {
namespace modelfit {

UnitSystem::UnitSystem(geom::SpherePoint const& position,
                       std::shared_ptr<const lsst::afw::image::PhotoCalib> photoCalib_, double flux) {
    auto scale = 1.0 * lsst::geom::arcseconds;
    auto cdMatrix = afw::geom::makeCdMatrix(scale);
    wcs = afw::geom::makeSkyWcs(geom::Point2D(0.0, 0.0), position, cdMatrix);
    photoCalib_ = validatePhotoCalib(photoCalib_);
    Scalar mag = photoCalib_->instFluxToMagnitude(flux);
    photoCalib = afw::image::makePhotoCalibFromCalibZeroPoint(std::pow(10.0, mag / 2.5), 0.0);
}

UnitSystem::UnitSystem(geom::SpherePoint const& position, Scalar mag) {
    auto scale = 1.0 * lsst::geom::arcseconds;
    auto cdMatrix = afw::geom::makeCdMatrix(scale);
    wcs = afw::geom::makeSkyWcs(geom::Point2D(0.0, 0.0), position, cdMatrix);
    photoCalib = afw::image::makePhotoCalibFromCalibZeroPoint(std::pow(10.0, mag / 2.5), 0.0);
}

std::shared_ptr<const lsst::afw::image::PhotoCalib> UnitSystem::validatePhotoCalib(
        std::shared_ptr<const lsst::afw::image::PhotoCalib> photoCalib_) {
    if (photoCalib_->getCalibrationMean() == 0.0) {
        return getDefaultPhotoCalib();
    } else {
        return photoCalib_;
    }
}

std::shared_ptr<const lsst::afw::image::PhotoCalib> UnitSystem::getDefaultPhotoCalib() {
    // Create a photoCalib object with a zero-point set to magnitude 27
    return afw::image::makePhotoCalibFromCalibZeroPoint(std::pow(10.0, 27.0 / 2.5), 0.0);
}

LocalUnitTransform::LocalUnitTransform(geom::Point2D const& sourcePixel, UnitSystem const& source,
                                       UnitSystem const& destination)
        : geometric(afw::geom::linearizeTransform(
                  *afw::geom::makeWcsPairTransform(*source.wcs, *destination.wcs), sourcePixel)),
          flux(destination.photoCalib->getInstFluxAtZeroMagnitude() /
               source.photoCalib->getInstFluxAtZeroMagnitude()),
          sb(flux / geometric.getLinear().computeDeterminant()) {}

}  // namespace modelfit
}  // namespace meas
}  // namespace lsst
