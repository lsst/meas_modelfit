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

#include "lsst/pex/config/python.h"

#include "lsst/meas/modelfit/DoubleShapeletPsfApprox.h"
#include "lsst/meas/modelfit/GeneralPsfFitter.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace modelfit {
namespace {

void declareDoubleShapelet(py::module &mod) {
    using Control = DoubleShapeletPsfApproxControl;
    using Algorithm = DoubleShapeletPsfApproxAlgorithm;

    using PyControl = py::class_<Control, std::shared_ptr<Control>>;
    using PyAlgorithm = py::class_<Algorithm, std::shared_ptr<Algorithm>, meas::base::SimpleAlgorithm>;

    PyControl clsControl(mod, "DoubleShapeletPsfApproxControl");
    clsControl.def(py::init<>());
    LSST_DECLARE_CONTROL_FIELD(clsControl, Control, innerOrder);
    LSST_DECLARE_CONTROL_FIELD(clsControl, Control, outerOrder);
    LSST_DECLARE_CONTROL_FIELD(clsControl, Control, radiusRatio);
    LSST_DECLARE_CONTROL_FIELD(clsControl, Control, peakRatio);
    LSST_DECLARE_CONTROL_FIELD(clsControl, Control, minRadius);
    LSST_DECLARE_CONTROL_FIELD(clsControl, Control, minRadiusDiff);
    LSST_DECLARE_CONTROL_FIELD(clsControl, Control, maxRadiusBoxFraction);
    LSST_DECLARE_NESTED_CONTROL_FIELD(clsControl, Control, optimizer);

    PyAlgorithm clsAlgorithm(mod, "DoubleShapeletPsfApproxAlgorithm");
    // wrap anonymous enum values as ints because we'll need to use them as ints
    clsAlgorithm.attr("Control") = clsControl;
    clsAlgorithm.attr("FAILURE") = py::cast(Algorithm::FAILURE);
    clsAlgorithm.attr("INVALID_POINT_FOR_PSF") = py::cast(Algorithm::INVALID_POINT_FOR_PSF);
    clsAlgorithm.attr("INVALID_MOMENTS") = py::cast(Algorithm::INVALID_MOMENTS);
    clsAlgorithm.attr("MAX_ITERATIONS") = py::cast(Algorithm::MAX_ITERATIONS);

    clsAlgorithm.def(py::init<Control const &, std::string const &, afw::table::Schema &>(), "ctrl"_a,
                     "name"_a, "schema"_a);
    clsAlgorithm.def_static("initializeResult", &Algorithm::initializeResult, "ctrl"_a);
    clsAlgorithm.def_static("fitMoments", &Algorithm::fitMoments, "result"_a, "ctrl"_a, "psfImage"_a);
    clsAlgorithm.def_static("makeObjective", &Algorithm::makeObjective, "moments"_a, "ctrl"_a, "psfImage"_a);
    clsAlgorithm.def_static("fitProfile", &Algorithm::fitProfile, "result"_a, "ctrl"_a, "psfImage"_a);
    clsAlgorithm.def_static("fitShapelets", &Algorithm::fitShapelets, "result"_a, "ctrl"_a, "psfImage"_a);
    clsAlgorithm.def("measure", &Algorithm::measure, "measRecord"_a, "exposure"_a);
    clsAlgorithm.def("fail", &Algorithm::fail, "measRecord"_a, "error"_a = nullptr);
}

void declareGeneral(py::module &mod) {
    using ComponentControl = GeneralPsfFitterComponentControl;
    using Control = GeneralPsfFitterControl;
    using Fitter = GeneralPsfFitter;
    using Algorithm = GeneralPsfFitterAlgorithm;

    using PyComponentControl = py::class_<ComponentControl, std::shared_ptr<ComponentControl>>;
    using PyControl = py::class_<Control, std::shared_ptr<Control>>;
    using PyFitter = py::class_<Fitter, std::shared_ptr<Fitter>>;
    using PyAlgorithm = py::class_<Algorithm, std::shared_ptr<Algorithm>, Fitter>;

    PyComponentControl clsComponentControl(mod, "GeneralPsfFitterComponentControl");
    clsComponentControl.def(py::init<int, double>(), "order"_a = 0, "radius"_a = 1.0);
    LSST_DECLARE_CONTROL_FIELD(clsComponentControl, ComponentControl, order);
    LSST_DECLARE_CONTROL_FIELD(clsComponentControl, ComponentControl, positionPriorSigma);
    LSST_DECLARE_CONTROL_FIELD(clsComponentControl, ComponentControl, ellipticityPriorSigma);
    LSST_DECLARE_CONTROL_FIELD(clsComponentControl, ComponentControl, radiusFactor);
    LSST_DECLARE_CONTROL_FIELD(clsComponentControl, ComponentControl, radiusPriorSigma);

    PyControl clsControl(mod, "GeneralPsfFitterControl");
    clsControl.def(py::init<>());
    LSST_DECLARE_NESTED_CONTROL_FIELD(clsControl, Control, inner);
    LSST_DECLARE_NESTED_CONTROL_FIELD(clsControl, Control, primary);
    LSST_DECLARE_NESTED_CONTROL_FIELD(clsControl, Control, wings);
    LSST_DECLARE_NESTED_CONTROL_FIELD(clsControl, Control, outer);
    LSST_DECLARE_NESTED_CONTROL_FIELD(clsControl, Control, optimizer);
    LSST_DECLARE_CONTROL_FIELD(clsControl, Control, defaultNoiseSigma);

    PyFitter clsFitter(mod, "GeneralPsfFitter");
    clsFitter.def(py::init<Control const &>(), "ctrl"_a);
    clsFitter.def("addFields", &Fitter::addFields, "schema"_a, "prefix"_a);
    clsFitter.def("getModel", &Fitter::getModel);
    clsFitter.def("getPrior", &Fitter::getPrior);
    clsFitter.def("adapt", &Fitter::adapt, "previousFit"_a, "previousModel"_a);
    // We use lambdas here because the C++ signature has an optional int*
    // argument not relevant for Python (which confuses pybind11).
    clsFitter.def("apply", [](Fitter const &self, afw::image::Image<Pixel> const &image,
                              afw::geom::ellipses::Quadrupole const &moments,
                              Scalar noiseSigma) { return self.apply(image, moments, noiseSigma); },
                  "image"_a, "moments"_a, "noiseSigma"_a = -1);
    clsFitter.def("apply", [](Fitter const &self, afw::image::Image<double> const &image,
                              afw::geom::ellipses::Quadrupole const &moments,
                              Scalar noiseSigma) { return self.apply(image, moments, noiseSigma); },
                  "image"_a, "moments"_a, "noiseSigma"_a = -1);
    clsFitter.def("apply", [](Fitter const &self, afw::image::Image<Pixel> const &image,
                              shapelet::MultiShapeletFunction const &initial,
                              Scalar noiseSigma) { return self.apply(image, initial, noiseSigma); },
                  "image"_a, "initial"_a, "noiseSigma"_a = -1);
    clsFitter.def("apply", [](Fitter const &self, afw::image::Image<double> const &image,
                              shapelet::MultiShapeletFunction const &initial,
                              Scalar noiseSigma) { return self.apply(image, initial, noiseSigma); },
                  "image"_a, "initial"_a, "noiseSigma"_a = -1);

    PyAlgorithm clsAlgorithm(mod, "GeneralPsfFitterAlgorithm");
    clsAlgorithm.attr("Control") = clsControl;
    clsAlgorithm.attr("FAILURE") = py::cast(Algorithm::FAILURE);
    clsAlgorithm.attr("MAX_INNER_ITERATIONS") = py::cast(Algorithm::MAX_INNER_ITERATIONS);
    clsAlgorithm.attr("MAX_OUTER_ITERATIONS") = py::cast(Algorithm::MAX_OUTER_ITERATIONS);
    clsAlgorithm.attr("EXCEPTION") = py::cast(Algorithm::EXCEPTION);
    clsAlgorithm.attr("CONTAINS_NAN") = py::cast(Algorithm::CONTAINS_NAN);
    clsAlgorithm.def(py::init<Control const &, afw::table::Schema &, std::string const &>(), "ctrl"_a,
                     "schema"_a, "prefix"_a);
    clsAlgorithm.def("getKey", &Algorithm::getKey);
    clsAlgorithm.def("measure",
                     (void (Algorithm::*)(afw::table::SourceRecord &, afw::image::Image<double> const &,
                                          shapelet::MultiShapeletFunction const &) const) &
                             Algorithm::measure,
                     "measRecord"_a, "image"_a, "initial"_a);
    clsAlgorithm.def("measure",
                     (void (Algorithm::*)(afw::table::SourceRecord &, afw::image::Image<double> const &,
                                          afw::geom::ellipses::Quadrupole const &) const) &
                             Algorithm::measure,
                     "measRecord"_a, "image"_a, "moments"_a);
    clsAlgorithm.def("fail", &Algorithm::fail, "measRecord"_a, "error"_a = nullptr);

    // MultiShapeletPsfLikelihood intentionally not exposed to Python.
}

PYBIND11_MODULE(psf, mod) {
    py::module::import("lsst.afw.image");
    py::module::import("lsst.afw.geom.ellipses");
    py::module::import("lsst.meas.base");
    py::module::import("lsst.shapelet");
    py::module::import("lsst.meas.modelfit.optimizer");

    declareDoubleShapelet(mod);
    declareGeneral(mod);
}

}
}
}
}  // namespace lsst::meas::modelfit::anonymous
