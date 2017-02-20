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

#include "lsst/meas/modelfit/MultiModel.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace modelfit {
namespace {

using PyMultiModel = py::class_<MultiModel, std::shared_ptr<MultiModel>, Model>;

PYBIND11_PLUGIN(multiModel) {

    py::module::import("lsst.meas.modelfit.model");

    py::module mod("multiModel");

    PyMultiModel cls(mod, "MultiModel");
    cls.def(
        py::init<ModelVector, MultiModel::NameVector const &>(),
        "components"_a, "prefixes"_a
    );
    cls.def("getComponents", &MultiModel::getComponents);

    // All other MultiModel methods are virtuals already inherited from
    // wrappers for Model.

    return mod.ptr();
}

}}}} // namespace lsst::meas::modelfit::anonymous
