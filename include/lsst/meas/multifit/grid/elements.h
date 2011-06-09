// -*- LSST-C++ -*-
/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010, 2011 LSST Corporation.
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

#ifndef LSST_MEAS_MULTIFIT_GRID_parameters
#define LSST_MEAS_MULTIFIT_GRID_parameters

#include "lsst/meas/multifit/definition/elements.h"

#include <boost/iterator/indirect_iterator.hpp>
#include <vector>
#include <set>
#include <Eigen/Core>

namespace lsst { namespace meas { namespace multifit { namespace grid {

template <SharedElementType E>
class SharedElement : public detail::SharedElementBase<E>, private boost::noncopyable {
public:
    
    // No ConstPtr typedef to make it clear that this class is strictly immutable.
    typedef boost::shared_ptr< SharedElement<E> > Ptr;
    typedef typename detail::SharedElementTraits<E>::Value Value;

    int const offset;

private:

    friend class Initializer;

    SharedElement(definition::SharedElement<E> const & definition, int offset_) : 
        detail::SharedElementBase<E>(definition), offset(offset_) {}

};


typedef SharedElement<POSITION> PositionElement;
typedef SharedElement<RADIUS> RadiusElement;
typedef SharedElement<ELLIPTICITY> EllipticityElement;

#ifndef SWIG

template <SharedElementType E>
class ElementArray {
    typedef boost::shared_ptr< SharedElement<T> > Ptr;
    typedef std::vector<Ptr> PtrVec;
    typedef typename PtrVec::const_iterator PtrIter;
public:

    typedef SharedElement<T> value_type;
    typedef Ptr pointer;
    typedef value_type const & reference;
    typedef reference const_reference;
    typedef std::ptrdiff_t difference_type;
    typedef std::size_t size_type;
    typedef boost::indirect_iterator<PtrIter> iterator;
    typedef iterator const_iterator;

    ElementArray() : _ptrVec() {}

    const_iterator begin() const { return _ptrVec.begin(); }
    const_iterator end() const { return _ptrVec.end(); }

    size_type size() const { return _ptrVec.size(); }

    bool empty() const { return _ptrVec.empty(); }

private:

    friend class Initializer;

    PtrVec _ptrVec;

};

template <SharedElementType E>
std::ostream & operator<<(std::ostream & os, SharedElement<E> const & component);

#endif

}}}} // namespace lsst::meas::multifit::grid

#endif // !LSST_MEAS_MULTIFIT_GRID_parameters
