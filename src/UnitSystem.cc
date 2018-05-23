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

namespace lsst { namespace meas { namespace modelfit {

UnitSystem::UnitSystem(afw::geom::SpherePoint const & position,
                       std::shared_ptr<const lsst::afw::image::Calib> calibIn,
                       double flux) {
    auto scale = 1.0*lsst::afw::geom::arcseconds;
    auto cdMatrix = afw::geom::makeCdMatrix(scale);
    wcs = afw::geom::makeSkyWcs(afw::geom::Point2D(0.0, 0.0), position, cdMatrix);
    calibIn = validateCalib(calibIn);
    Scalar mag = calibIn->getMagnitude(flux);
    PTR(afw::image::Calib) calib_ = std::make_shared<afw::image::Calib>();
    calib_->setFluxMag0(std::pow(10.0, mag/2.5));
    calib = calib_;
}

UnitSystem::UnitSystem(afw::geom::SpherePoint const & position, Scalar mag) {
    auto scale = 1.0*lsst::afw::geom::arcseconds;
    auto cdMatrix = afw::geom::makeCdMatrix(scale);
    wcs = afw::geom::makeSkyWcs(afw::geom::Point2D(0.0, 0.0), position, cdMatrix);
    PTR(afw::image::Calib) calib_ = std::make_shared<afw::image::Calib>();
    calib_->setFluxMag0(std::pow(10.0, mag/2.5));
    calib = calib_;
}

std::shared_ptr<const lsst::afw::image::Calib> UnitSystem::validateCalib(
    std::shared_ptr<const lsst::afw::image::Calib> calib_) {
    if (calib_->getFluxMag0().first == 0.0){
        return getDefaultCalib();
    } else{
        return calib_;
    }
}

std::shared_ptr<const lsst::afw::image::Calib> UnitSystem::getDefaultCalib() {
    // Create a default calib object with a zero-point set to magnitude 27
    static std::shared_ptr<const lsst::afw::image::Calib> tmp =
                                                        std::make_shared<const lsst::afw::image::Calib>(
                                                        std::pow(10.0, 27.0/2.5)
                                                        );
    return tmp;
}

LocalUnitTransform::LocalUnitTransform(
    afw::geom::Point2D const & sourcePixel,
    UnitSystem const & source,
    UnitSystem const & destination
) :
    geometric(afw::geom::linearizeTransform(*afw::geom::makeWcsPairTransform(*source.wcs, *destination.wcs),
                                            sourcePixel)),
    flux(destination.calib->getFluxMag0().first / source.calib->getFluxMag0().first),
    sb(flux / geometric.getLinear().computeDeterminant())
{}

}}} // namespace lsst::meas::modelfit
