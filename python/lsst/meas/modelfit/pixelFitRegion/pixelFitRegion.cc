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

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "lsst/pex/config/pybind11.h"

#include "lsst/meas/modelfit/PixelFitRegion.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace modelfit {
namespace {

using PyPixelFitRegionControl = py::class_<PixelFitRegionControl, std::shared_ptr<PixelFitRegionControl>>;
using PyPixelFitRegion = py::class_<PixelFitRegion, std::shared_ptr<PixelFitRegion>>;

PYBIND11_PLUGIN(pixelFitRegion) {

    py::module::import("lsst.afw.image");
    py::module::import("lsst.afw.detection");
    py::module::import("lsst.afw.geom.ellipses");

    py::module mod("pixelFitRegion");

    using Control = PixelFitRegionControl;

    PyPixelFitRegionControl clsControl(mod, "PixelFitRegionControl");
    clsControl.def(py::init<>());
    LSST_DECLARE_CONTROL_FIELD(clsControl, Control, nKronRadii);
    LSST_DECLARE_CONTROL_FIELD(clsControl, Control, nPsfSigmaMin);
    LSST_DECLARE_CONTROL_FIELD(clsControl, Control, nPsfSigmaGrow);
    LSST_DECLARE_CONTROL_FIELD(clsControl, Control, nFitRadiiMin);
    LSST_DECLARE_CONTROL_FIELD(clsControl, Control, nFitRadiiMax);
    LSST_DECLARE_CONTROL_FIELD(clsControl, Control, maxArea);
    LSST_DECLARE_CONTROL_FIELD(clsControl, Control, badMaskPlanes);
    LSST_DECLARE_CONTROL_FIELD(clsControl, Control, maxBadPixelFraction);

    PyPixelFitRegion cls(mod, "PixelFitRegion");
    cls.def(
        py::init<
            Control const &,
            afw::geom::ellipses::Quadrupole const &,
            afw::geom::ellipses::Quadrupole const &,
            Scalar,
            int
        >(),
        "ctrl"_a, "moments"_a, "psfMoments"_a, "kronRadius"_a, "footprintArea"_a
    );
    cls.def(
        py::init<Control const &, afw::geom::ellipses::Quadrupole const &>(),
        "ctrl"_a, "ellipse"_a
    );
    cls.def("applyEllipse", &PixelFitRegion::applyEllipse, "deconvolved"_a, "psfMoments"_a);
    cls.def("applyMask", &PixelFitRegion::applyMask, "mask"_a, "center"_a);
    // Data members are intentionally read-only from the Python side;
    // they should only be set by the constructor and apply methods.
    cls.def_readonly("ellipse", &PixelFitRegion::ellipse);
    cls.def_readonly("footprint", &PixelFitRegion::footprint);
    cls.def_readonly("usedFootprintArea", &PixelFitRegion::usedFootprintArea);
    cls.def_readonly("usedPsfArea", &PixelFitRegion::usedPsfArea);
    cls.def_readonly("maxArea", &PixelFitRegion::maxArea);
    cls.def_readonly("maxBadPixelFraction", &PixelFitRegion::maxBadPixelFraction);
    cls.def_readonly("usedMinEllipse", &PixelFitRegion::usedMinEllipse);
    cls.def_readonly("usedMaxEllipse", &PixelFitRegion::usedMaxEllipse);

    return mod.ptr();
}

}}}} // namespace lsst::meas::modelfit::anonymous
