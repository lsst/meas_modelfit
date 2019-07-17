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

#ifndef LSST_MEAS_MODELFIT_UnitSystem_h_INCLUDED
#define LSST_MEAS_MODELFIT_UnitSystem_h_INCLUDED

#include "lsst/afw/image/Exposure.h"
#include "lsst/geom/SpherePoint.h"
#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/afw/image/PhotoCalib.h"
#include "lsst/geom/AffineTransform.h"

#include "lsst/meas/modelfit/common.h"

namespace lsst {
namespace meas {
namespace modelfit {

/**
 *  @brief A simple struct that combines a Wcs and a PhotoCalib.
 */
struct UnitSystem {
    std::shared_ptr<afw::geom::SkyWcs const> wcs;
    std::shared_ptr<afw::image::PhotoCalib const> photoCalib;

    /**
     *  @brief Construct a "standard" UnitSystem
     *
     *  This constructs a UnitSystem with a TAN Wcs centered on the given position, with flux units
     *  set such that unit flux is the given magnitude.  See @ref modelfitUnits for an explanation
     *  of why we frequently use this system.
     */
    UnitSystem(geom::SpherePoint const& position,
               std::shared_ptr<const afw::image::PhotoCalib> photoCalib, double flux);
    UnitSystem(geom::SpherePoint const& position, Scalar mag);

    /// Construct a UnitSystem from a given Wcs and PhotoCalib
    UnitSystem(std::shared_ptr<afw::geom::SkyWcs const> wcs_,
               std::shared_ptr<afw::image::PhotoCalib const> photoCalib)
            : wcs(wcs_), photoCalib(validatePhotoCalib(photoCalib)) {}

    /// Construct a UnitSystem by extracting the Wcs and PhotoCalib from an Exposure (implicit)
    template <typename T>
    UnitSystem(afw::image::Exposure<T> const& exposure)
            : wcs(exposure.getWcs()), photoCalib(validatePhotoCalib(exposure.getCalib())) {}

private:
    std::shared_ptr<const afw::image::PhotoCalib> validatePhotoCalib(
            std::shared_ptr<const afw::image::PhotoCalib> photoCalib);
    static std::shared_ptr<const afw::image::PhotoCalib> getDefaultPhotoCalib();
};

/**
 *  @brief A local mapping between two UnitSystems
 *
 *  LocalUnitTransform is "local" because it linearizes the Wcs and evaluates the PhotoCalib transform
 *  at a particular predifined point, allowing it to represent the geometric transform as an
 *  AffineTransform and the photometric transform as a simple scaling.
 */
struct LocalUnitTransform {
    /// Maps source pixel coordinates to destination pixel coordinates
    geom::AffineTransform geometric;

    /// Multiply source fluxes by this to get destination fluxes
    double flux;

    /// Multiply source surface brightnesses by this to get destination surface brightnesses
    double sb;

    LocalUnitTransform(geom::Point2D const& sourcePixel, UnitSystem const& source,
                       UnitSystem const& destination);

    /// Construct an identity transform for both geometry and flux.
    LocalUnitTransform() : geometric(), flux(1.0), sb(1.0) {}
};

}  // namespace modelfit
}  // namespace meas
}  // namespace lsst

#endif  // !LSST_MEAS_MODELFIT_UnitSystem_h_INCLUDED
