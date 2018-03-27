/*
 * This file is part of meas_modelfit.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
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
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "numpy/arrayobject.h"
#include "ndarray/pybind11.h"
#include "ndarray/eigen.h"

#include "lsst/meas/modelfit/RegularizedMoments.h"

namespace py = pybind11;
using namespace py::literals;

namespace lsst {
namespace meas {
namespace modelfit {

namespace {
    using PyMomentsModel = py::class_<MomentsModel, std::shared_ptr<MomentsModel>>;

PYBIND11_PLUGIN(regularizedMoments) {
    if (_import_array() < 0) {
        PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import");
        return nullptr;
    }

    py::module mod("regularizedMoments");
    PyMomentsModel cls(mod, "MomentsModel");

    cls.def(py::init<MomentsModel::Moments const &>(), "weights"_a);

    cls.def("at", &MomentsModel::at);
    cls.def("computeValues", &MomentsModel::computeValues);
    cls.def("computeJacobian", &MomentsModel::computeJacobian);

    // Wrap the test functions for the anonymous classes
    mod.def("testNorm", &testNorm, "tol"_a=1e-6);
    mod.def("testAlphaX", &testAlphaX, "tol"_a=1e-6);
    mod.def("testAlphaY", &testAlphaY, "tol"_a=1e-6);
    mod.def("testBetaX", &testBetaX, "tol"_a=1e-6);
    mod.def("testBetaXY", &testBetaXY, "tol"_a=1e-6);
    mod.def("testBetaY", &testBetaY, "tol"_a=1e-6);

    return mod.ptr();
}
} // end anonymous
}}} // end lsst::meas::modelfit
