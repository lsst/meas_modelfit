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

#ifndef LSST_MEAS_MULTIFIT_UnitSystem_h_INCLUDED
#define LSST_MEAS_MULTIFIT_UnitSystem_h_INCLUDED

#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/Calib.h"
#include "lsst/afw/geom/AffineTransform.h"
#include "lsst/afw/geom/XYTransform.h"

#include "lsst/meas/multifit/common.h"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief A simple struct that combines a Wcs and a Calib.
 */
struct UnitSystem {
    PTR(afw::image::Wcs const) wcs;
    PTR(afw::image::Calib const) calib;

    /**
     *  @brief Construct a "standard" UnitSystem
     *
     *  This constructs a UnitSystem with a TAN Wcs centered on the given position, with flux units
     *  set such that unit flux is the given magnitude.  See @ref multifitUnits for an explanation
     *  of why we frequently use this system.
     */
    UnitSystem(afw::coord::Coord const & position, Scalar mag);

    /// Construct a UnitSystem from a give Wcs and Calib
    UnitSystem(PTR(afw::image::Wcs const) wcs_, PTR(afw::image::Calib const) calib_) :
        wcs(wcs_), calib(calib_)
    {}

    /// Construct a UnitSystem by extracting the Wcs and Calib from an Exposure (implicit)
    template <typename T>
    UnitSystem(afw::image::Exposure<T> const & exposure) :
        wcs(exposure.getWcs()), calib(exposure.getCalib())
    {}
};

/**
 *  @brief A local mapping between two UnitSystems
 *
 *  LocalUnitTransform is "local" because it linearizes the Wcs and evaluates the Calib transform
 *  at a particular predifined point, allowing it to represent the geometric transform as an
 *  AffineTransform and the photometric transform as a simple scaling.
 */
struct LocalUnitTransform {

    /// Maps source pixel coordinates to destination pixel coordinates
    afw::geom::AffineTransform geometric;

    /// Multiply source fluxes by this to get destination fluxes
    double flux;

    /// Multiply source surface brightnesses by this to get destination surface brightnesses
    double sb;

    LocalUnitTransform(
        afw::coord::Coord const & position,
        UnitSystem const & source,
        UnitSystem const & destination
    );

};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_UnitSystem_h_INCLUDED
