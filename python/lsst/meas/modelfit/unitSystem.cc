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

#include "lsst/meas/modelfit/UnitSystem.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace modelfit {
namespace {

using PyUnitSystem = py::class_<UnitSystem, std::shared_ptr<UnitSystem>>;
using PyLocalUnitTransform = py::class_<LocalUnitTransform, std::shared_ptr<LocalUnitTransform>>;

PYBIND11_PLUGIN(unitSystem) {
    py::module::import("lsst.afw.coord");
    py::module::import("lsst.afw.image");

    py::module mod("unitSystem");

    // Data member wrappers in this module are intentionally read-only;
    // setting them in Python should never be necessary and hence always
    // represents a mistake we want to fail on as early as possible.

    PyUnitSystem clsUnitSystem(mod, "UnitSystem");
    clsUnitSystem.def_readonly("wcs", &UnitSystem::wcs);
    clsUnitSystem.def_readonly("calib", &UnitSystem::calib);
    clsUnitSystem.def(py::init<afw::coord::Coord const &, std::shared_ptr<afw::image::Calib const>, double>(),
                      "position"_a, "calibIn"_a, "flux"_a);
    clsUnitSystem.def(py::init<afw::coord::Coord const &, Scalar>(), "position"_a, "mag"_a);
    clsUnitSystem.def(
            py::init<std::shared_ptr<afw::image::Wcs const>, std::shared_ptr<afw::image::Calib const>>(),
            "wcs"_a, "calib"_a);
    clsUnitSystem.def(py::init<afw::image::Exposure<float> const &>(), "exposure"_a);
    clsUnitSystem.def(py::init<afw::image::Exposure<double> const &>(), "exposure"_a);

    PyLocalUnitTransform clsLocalUnitTransform(mod, "LocalUnitTransform");
    clsLocalUnitTransform.def_readonly("geometric", &LocalUnitTransform::geometric);
    clsLocalUnitTransform.def_readonly("flux", &LocalUnitTransform::flux);
    clsLocalUnitTransform.def_readonly("sb", &LocalUnitTransform::sb);
    clsLocalUnitTransform.def(py::init<afw::coord::Coord const &, UnitSystem const &, UnitSystem const &>(),
                              "position"_a, "source"_a, "destination"_a);
    clsLocalUnitTransform.def(py::init<>());

    return mod.ptr();
}
}
}
}
}  // namespace lsst::meas::modelfit::anonymous
