// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2015 LSST/AURA.
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

#ifndef LSST_MEAS_MODELFIT_DebugDump_h_INCLUDED
#define LSST_MEAS_MODELFIT_DebugDump_h_INCLUDED

#include "lsst/base.h"
#include "ndarray_fwd.h"

namespace lsst { namespace meas { namespace modelfit {

class DebugDump {
public:

    explicit DebugDump(std::string const & name) : _dict(nullptr) {
#ifdef LSST_DEBUG_DUMP
        _initialize(name);
#endif
    }

    template <typename T>
    void operator()(std::string const & name, T const & value) {
#ifdef LSST_DEBUG_DUMP
        _dump(name, value);
#endif
    }

    template <typename T, int N, int C>
    void copy(std::string const & name, ndarray::Array<T const,N,C> const & value) {
#ifdef LSST_DEBUG_DUMP
        _dump(name, ndarray::copy(value).shallow());
#endif
    }

    template <typename T, int N, int C>
    void copy(std::string const & name, ndarray::Array<T,N,C> const & value) {
#ifdef LSST_DEBUG_DUMP
        _dump(name, ndarray::copy(value).shallow());
#endif
    }

    DebugDump(DebugDump const &) = delete;
    DebugDump(DebugDump &&) = delete;
    DebugDump & operator=(DebugDump const &) = delete;
    DebugDump & operator=(DebugDump &&) = delete;

    ~DebugDump() {
#ifdef LSST_DEBUG_DUMP
        _destroy();
#endif
    }

private:

    class Impl;

    void _initialize(std::string const & name);

    void _destroy();

    template <typename T>
    void _dump(std::string const & name, T const & scalar);

    void * _dict;
};

}}} // namespace lsst::meas:modelfit

#endif // !LSST_MEAS_MODELFIT_DebugDump_h_INCLUDED
