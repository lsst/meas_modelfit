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

%{

#include "ndarray/swig/numpy.h"
#include "lsst/meas/modelfit/DebugDump.h"

namespace lsst { namespace meas { namespace modelfit {

void DebugDump::_initialize(std::string const & name) {
    _dict = PyDict_New();
    PyObject * module = PyImport_ImportModule("lsst.meas.modelfit");
    if (!module) {
        PyErr_Warn(nullptr, "Failed to import module for debug dumps.");
    }
    PyObject * registry = PyObject_GetAttrString(module, "debugRegistry");
    if (!registry) {
        PyErr_Clear();
        registry = PyDict_New();
        Py_INCREF(registry);
        PyObject_SetAttrString(module, "debugRegistry", registry);
    }
    PyObject * dict = reinterpret_cast<PyObject*>(_dict);
    Py_INCREF(dict);
    PyDict_SetItemString(registry, name.c_str(), dict);
    Py_DECREF(registry);
    Py_DECREF(module);
}

void DebugDump::_destroy() {
    if (_dict) {
        PyObject * p = reinterpret_cast<PyObject *>(_dict);
        Py_DECREF(p);
    }
}

template <typename T>
void DebugDump::_dump(std::string const & name, T const & value) {
    PyObject * result = ndarray::PyConverter<T>::toPython(value);
    if (!result) {
        PyErr_Warn(nullptr, "Failed to dump debug object: could not convert to Python.");
        return;
    }
    PyDict_SetItemString(reinterpret_cast<PyObject*>(_dict), name.c_str(), result);
}

template void DebugDump::_dump(std::string const & name, int const &);
template void DebugDump::_dump(std::string const & name, float const &);
template void DebugDump::_dump(std::string const & name, double const &);

template void DebugDump::_dump(std::string const & name, ndarray::Array<float,1,0> const &);
template void DebugDump::_dump(std::string const & name, ndarray::Array<float,1,1> const &);
template void DebugDump::_dump(std::string const & name, ndarray::Array<float,2,-1> const &);
template void DebugDump::_dump(std::string const & name, ndarray::Array<float,2,-2> const &);
template void DebugDump::_dump(std::string const & name, ndarray::Array<float,2,0> const &);
template void DebugDump::_dump(std::string const & name, ndarray::Array<float,2,1> const &);
template void DebugDump::_dump(std::string const & name, ndarray::Array<float,2,2> const &);

template void DebugDump::_dump(std::string const & name, ndarray::Array<float const,1,0> const &);
template void DebugDump::_dump(std::string const & name, ndarray::Array<float const,1,1> const &);
template void DebugDump::_dump(std::string const & name, ndarray::Array<float const,2,-1> const &);
template void DebugDump::_dump(std::string const & name, ndarray::Array<float const,2,-2> const &);
template void DebugDump::_dump(std::string const & name, ndarray::Array<float const,2,0> const &);
template void DebugDump::_dump(std::string const & name, ndarray::Array<float const,2,1> const &);
template void DebugDump::_dump(std::string const & name, ndarray::Array<float const,2,2> const &);

template void DebugDump::_dump(std::string const & name, ndarray::Array<double,1,0> const &);
template void DebugDump::_dump(std::string const & name, ndarray::Array<double,1,1> const &);
template void DebugDump::_dump(std::string const & name, ndarray::Array<double,2,-1> const &);
template void DebugDump::_dump(std::string const & name, ndarray::Array<double,2,-2> const &);
template void DebugDump::_dump(std::string const & name, ndarray::Array<double,2,0> const &);
template void DebugDump::_dump(std::string const & name, ndarray::Array<double,2,1> const &);
template void DebugDump::_dump(std::string const & name, ndarray::Array<double,2,2> const &);

template void DebugDump::_dump(std::string const & name, ndarray::Array<double const,1,0> const &);
template void DebugDump::_dump(std::string const & name, ndarray::Array<double const,1,1> const &);
template void DebugDump::_dump(std::string const & name, ndarray::Array<double const,2,-1> const &);
template void DebugDump::_dump(std::string const & name, ndarray::Array<double const,2,-2> const &);
template void DebugDump::_dump(std::string const & name, ndarray::Array<double const,2,0> const &);
template void DebugDump::_dump(std::string const & name, ndarray::Array<double const,2,1> const &);
template void DebugDump::_dump(std::string const & name, ndarray::Array<double const,2,2> const &);

template void DebugDump::_dump(std::string const & name, Eigen::Matrix<float,Eigen::Dynamic,1> const &);
template void DebugDump::_dump(std::string const & name, Eigen::Matrix<float,1,Eigen::Dynamic> const &);
template void DebugDump::_dump(std::string const & name, Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> const &);
template void DebugDump::_dump(std::string const & name, Eigen::Matrix<double,Eigen::Dynamic,1> const &);
template void DebugDump::_dump(std::string const & name, Eigen::Matrix<double,1,Eigen::Dynamic> const &);
template void DebugDump::_dump(std::string const & name, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> const &);


}}} // namespace lsst::meas::modelfit

%}
