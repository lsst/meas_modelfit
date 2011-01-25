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

#ifndef LSST_MEAS_MULTIFIT_GRID_parameters
#define LSST_MEAS_MULTIFIT_GRID_parameters

#include "multifit/definition/parameters.hpp"

#include <boost/iterator/indirect_iterator.hpp>
#include <vector>
#include <set>
#include <Eigen/Core>

namespace lsst { namespace meas { namespace multifit {

class Grid;

namespace grid {

template <parameters::Enum E>
class ParameterComponent : public definition::ParameterComponent<E> {
public:

    typedef definition::ParameterComponent<E> Base;

    ParameterComponent(Base const & definition, int offset_) : Base(definition), offset(offset_) {}

    ParameterComponent(ParameterComponent const & other) : Base(other), offset(other.offset) {}

    int const offset;

    void writeParameters(double * parameters) const { this->_writeParameters(parameters + offset); }

    typename Base::Ptr makeDefinition() const { return this->copy(); }

    typename Base::Ptr makeDefinition(double const * parameters) {
        typename Base::Ptr result = this->copy();
        result->_readParameters(parameters + offset);
        return result;
    }

};

template <>
class ParameterComponent<parameters::RADIUS> : public definition::RadiusComponent {
public:

    typedef definition::RadiusComponent Base;

    ParameterComponent(Base const & definition, int offset_) : Base(definition), offset(offset_) {}

    ParameterComponent(ParameterComponent const & other) : Base(other), offset(other.offset) {}

    int const offset;

    std::set<definition::EllipticityComponent::Ptr> associated_ellipticities;

    void writeParameters(double * parameters) const { this->_writeParameters(parameters + offset); }

    Base::Ptr makeDefinition() const { return this->copy(); }

    Base::Ptr makeDefinition(double const * parameters) {
        Base::Ptr result = this->copy();
        result->_readParameters(parameters + offset);
        return result;
    }

};

typedef ParameterComponent<parameters::POSITION> PositionComponent;
typedef ParameterComponent<parameters::RADIUS> RadiusComponent;
typedef ParameterComponent<parameters::ELLIPTICITY> EllipticityComponent;

template <typename T>
class ComponentArray {
    typedef boost::shared_ptr<T> Ptr;
    typedef std::vector<Ptr> PtrVec;
    typedef typename PtrVec::const_iterator PtrIter;
public:

    typedef T value_type;
    typedef T * pointer;
    typedef T & reference;
    typedef T const & const_reference;
    typedef std::ptrdiff_t difference_type;
    typedef std::size_t size_type;
    typedef boost::indirect_iterator<PtrIter> iterator;
    typedef boost::indirect_iterator<PtrIter,T const> const_iterator;

    ComponentArray() : _ptr_vec() {}

    const_iterator begin() const { return _ptr_vec.begin(); }
    const_iterator end() const { return _ptr_vec.end(); }

    size_type size() const { return _ptr_vec.size(); }

    bool empty() const { return _ptr_vec.empty(); }

private:

    friend class multifit::Grid;

    PtrVec _ptr_vec;

};

}}}} // namespace lsst::meas::multifit::grid

#endif // !LSST_MEAS_MULTIFIT_GRID_parameters
