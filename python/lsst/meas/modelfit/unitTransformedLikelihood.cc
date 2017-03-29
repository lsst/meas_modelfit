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

#include "numpy/arrayobject.h"
#include "ndarray/pybind11.h"

#include "lsst/pex/config/python.h"
#include "lsst/meas/modelfit/UnitTransformedLikelihood.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace modelfit {
namespace {

using PyUnitTransformedLikelihoodControl =
        py::class_<UnitTransformedLikelihoodControl, std::shared_ptr<UnitTransformedLikelihoodControl>>;

using PyEpochFootprint = py::class_<EpochFootprint, std::shared_ptr<EpochFootprint>>;

using PyUnitTransformedLikelihood =
        py::class_<UnitTransformedLikelihood, std::shared_ptr<UnitTransformedLikelihood>, Likelihood>;

PYBIND11_PLUGIN(unitTransformedLikelihood) {
    py::module::import("lsst.afw.geom.ellipses");
    py::module::import("lsst.afw.detection");
    py::module::import("lsst.afw.image");
    py::module::import("lsst.meas.modelfit.model");
    py::module::import("lsst.meas.modelfit.likelihood");
    py::module::import("lsst.meas.modelfit.unitSystem");

    py::module mod("unitTransformedLikelihood");

    if (_import_array() < 0) {
        PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import");
        return nullptr;
    }

    PyUnitTransformedLikelihoodControl clsControl(mod, "UnitTransformedLikelihoodControl");
    LSST_DECLARE_CONTROL_FIELD(clsControl, UnitTransformedLikelihoodControl, usePixelWeights);
    LSST_DECLARE_CONTROL_FIELD(clsControl, UnitTransformedLikelihoodControl, weightsMultiplier);
    clsControl.def(py::init<bool>(), "usePixelWeights"_a = false);

    PyEpochFootprint clsEpochFootprint(mod, "EpochFootprint");
    clsEpochFootprint.def(py::init<afw::detection::Footprint const &, afw::image::Exposure<Pixel> const &,
                                   shapelet::MultiShapeletFunction const &>(),
                          "footprint"_a, "exposure"_a, "psf"_a);
    clsEpochFootprint.def_readonly("footprint", &EpochFootprint::footprint);
    clsEpochFootprint.def_readonly("exposure", &EpochFootprint::exposure);
    clsEpochFootprint.def_readonly("psf", &EpochFootprint::psf);

    PyUnitTransformedLikelihood clsUnitTransformedLikelihood(mod, "UnitTransformedLikelihood");
    clsUnitTransformedLikelihood.def(
            py::init<std::shared_ptr<Model>, ndarray::Array<Scalar const, 1, 1> const &, UnitSystem const &,
                     afw::coord::Coord const &, afw::image::Exposure<Pixel> const &,
                     afw::detection::Footprint const &, shapelet::MultiShapeletFunction const &,
                     UnitTransformedLikelihoodControl const &>(),
            "model"_a, "fixed"_a, "fitSys"_a, "position"_a, "exposure"_a, "footprint"_a, "psf"_a, "ctrl"_a);
    clsUnitTransformedLikelihood.def(
            py::init<std::shared_ptr<Model>, ndarray::Array<Scalar const, 1, 1> const &, UnitSystem const &,
                     afw::coord::Coord const &, std::vector<std::shared_ptr<EpochFootprint>> const &,
                     UnitTransformedLikelihoodControl const &>(),
            "model"_a, "fixed"_a, "fitSys"_a, "position"_a, "epochFootprintList"_a, "ctrl"_a);

    return mod.ptr();
}
}
}
}
}  // namespace lsst::meas::modelfit::anonymous
