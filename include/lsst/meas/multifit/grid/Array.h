// -*- LSST-C++ -*-
/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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

#ifndef LSST_MEAS_MULTIFIT_GRID_Array
#define LSST_MEAS_MULTIFIT_GRID_Array

#include <boost/noncopyable.hpp>
#include <boost/iterator/iterator_adaptor.hpp>

namespace lsst { namespace meas { namespace multifit { 

class Grid;

namespace grid {

template <typename T>
class Array : private boost::noncopyable {
public:

    typedef T value_type;
    typedef T * pointer;
    typedef T & reference;
    typedef T const & const_reference;
    typedef std::ptrdiff_t difference_type;
    typedef std::size_t size_type;
    typedef T * iterator;
    typedef T const * const_iterator;

    Array() : _first(0), _last(0) {}

    const_iterator begin() const { return _first; }
    const_iterator end() const { return _last; }

    size_type size() const { return _last - _first; }

    bool empty() const { return _last == _first; }

private:

    friend class multifit::Grid;

    explicit Array(pointer first, pointer last) :
        _first(first), _last(last) {}

    pointer _first;
    pointer _last;
};

template <typename T>
T const & find(Array<T> const & array, ID const id);

}}}} // namespace lsst::meas::multifit::grid

#endif // !LSST_MEAS_MULTIFIT_GRID_Array
