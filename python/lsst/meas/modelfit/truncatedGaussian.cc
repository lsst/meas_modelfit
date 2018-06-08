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

#include "ndarray/pybind11.h"

#include "lsst/meas/modelfit/TruncatedGaussian.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace modelfit {
namespace {

using Sampler = TruncatedGaussianSampler;
using Evaluator = TruncatedGaussianEvaluator;
using LogEvaluator = TruncatedGaussianLogEvaluator;

using PyTruncatedGaussian = py::class_<TruncatedGaussian, std::shared_ptr<TruncatedGaussian>>;
using PySampler = py::class_<Sampler, std::shared_ptr<Sampler>>;
using PyEvaluator = py::class_<Evaluator, std::shared_ptr<Evaluator>>;
using PyLogEvaluator = py::class_<LogEvaluator, std::shared_ptr<LogEvaluator>>;

// Shared wrapper code for the TruncatedGaussianLogEvaluator and
// TruncatedGaussianEvaluator classes, which have the exact same interface.
// 'name' should be one of ("Evaluator", "LogEvaluator").
template <typename Class, typename PyClass>
static PyClass declareEvaluator(py::module &mod, std::string const &name) {
    PyClass cls(mod, ("TruncatedGaussian" + name).c_str());
    cls.def(py::init<TruncatedGaussian const &>(), "parent"_a);
    cls.def("__call__",
            (Scalar (Class::*)(ndarray::Array<Scalar const, 1, 1> const &) const) & Class::operator(),
            "alpha"_a);
    cls.def("__call__", (void (Class::*)(ndarray::Array<Scalar const, 2, 1> const &,
                                         ndarray::Array<Scalar, 1, 1> const &) const) &
                                Class::operator(),
            "alpha"_a, "output"_a);
    // Third overload of operator() is just an Eigen version of the ndarray one, so it's
    // redundant in Python.
    return cls;
}

PYBIND11_PLUGIN(truncatedGaussian) {
    py::module::import("lsst.afw.math");

    py::module mod("truncatedGaussian");

    PyTruncatedGaussian cls(mod, "TruncatedGaussian");
    py::enum_<TruncatedGaussian::SampleStrategy>(cls, "SampleStrategy")
            .value("DIRECT_WITH_REJECTION", TruncatedGaussian::DIRECT_WITH_REJECTION)
            .value("ALIGN_AND_WEIGHT", TruncatedGaussian::ALIGN_AND_WEIGHT)
            .export_values();
    cls.def_static("fromSeriesParameters", &TruncatedGaussian::fromSeriesParameters, "q0"_a, "gradient"_a,
                   "hessian"_a);
    cls.def_static("fromStandardParameters", &TruncatedGaussian::fromStandardParameters, "mean"_a,
                   "covariance"_a);
    cls.def("sample", (Sampler (TruncatedGaussian::*)(TruncatedGaussian::SampleStrategy) const) &
                              TruncatedGaussian::sample,
            "strategy"_a);
    cls.def("sample", (Sampler (TruncatedGaussian::*)(Scalar) const) & TruncatedGaussian::sample,
            "minRejectionEfficiency"_a = 0.1);
    cls.def("evaluateLog", &TruncatedGaussian::evaluateLog);
    cls.def("evaluate", &TruncatedGaussian::evaluate);
    cls.def("getDim", &TruncatedGaussian::getDim);
    cls.def("maximize", &TruncatedGaussian::maximize);
    cls.def("getUntruncatedFraction", &TruncatedGaussian::getUntruncatedFraction);
    cls.def("getLogPeakAmplitude", &TruncatedGaussian::getLogPeakAmplitude);
    cls.def("getLogIntegral", &TruncatedGaussian::getLogIntegral);

    cls.attr("LogEvaluator") = declareEvaluator<LogEvaluator, PyLogEvaluator>(mod, "LogEvaluator");
    cls.attr("Evaluator") = declareEvaluator<Evaluator, PyEvaluator>(mod, "Evaluator");

    PySampler clsSampler(mod, "TruncatedGaussianSampler");
    clsSampler.def(py::init<TruncatedGaussian const &, TruncatedGaussian::SampleStrategy>(), "parent"_a,
                   "strategy"_a);
    clsSampler.def("__call__",
                   (Scalar (Sampler::*)(afw::math::Random &, ndarray::Array<Scalar, 1, 1> const &) const) &
                           Sampler::operator(),
                   "rng"_a, "alpha"_a);
    clsSampler.def("__call__", (void (Sampler::*)(afw::math::Random &, ndarray::Array<Scalar, 2, 1> const &,
                                                  ndarray::Array<Scalar, 1, 1> const &, bool) const) &
                                       Sampler::operator(),
                   "rng"_a, "alpha"_a, "weights"_a, "multiplyWeights"_a = false);

    cls.attr("Sampler") = clsSampler;

    return mod.ptr();
}
}
}
}
}  // namespace lsst::meas::modelfit::anonymous
