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

#ifndef LSST_MEAS_MULTIFIT_DEFINITION_parameters
#define LSST_MEAS_MULTIFIT_DEFINITION_parameters

#include "lsst/meas/multifit/constants.h"
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

namespace lsst { namespace meas { namespace multifit {

namespace grid {
template <ParameterType E> class ParameterComponent;
} // namespace grid

namespace definition {

/**
 *  @brief A base class for parameter management classes.
 *
 *  ParameterComponents can be shared by multiple objects, tying
 *  those parameters together (such as the bulge and disk of a galaxy,
 *  for instance).
 *
 *  A ParameterComponent is 'active' if its parameters are allowed
 *  to vary in an optimization procedure; it should maintain enough
 *  state to be able to compute its products (which are 
 *  subclass-dependent) even when inactive.
 */
template <typename Derived, typename Value_, int size_>
class ParameterComponentBase {
public:

    typedef boost::shared_ptr<Derived> Ptr;
    typedef Value_ Value;

    static int const SIZE = size_;

    bool active;

    Value const & getValue() const { return _value; }

    Ptr copy() const { return boost::make_shared<Derived>(static_cast<Derived const &>(*this)); }

protected:

    explicit ParameterComponentBase(Value const & value) : 
        active(true), _value(value) {}

    ParameterComponentBase(ParameterComponentBase const & other) :
        active(other.active), _value(other._value) {}

    Value _value;
};

template <ParameterType P> class ParameterComponent;

/**
 *  @brief Parameters representing a position.
 *
 *  A position is represented as an offset from a reference point, which is stored within the
 *  position.
 */
template <>
class ParameterComponent<POSITION>
    : public ParameterComponentBase<ParameterComponent<POSITION>,lsst::afw::geom::Extent2D,2>
{
    typedef ParameterComponentBase<ParameterComponent<POSITION>,lsst::afw::geom::Extent2D,2> Base;
public:

    explicit ParameterComponent(lsst::afw::geom::Point2D const & reference, Value const & value = Value()) :
        Base(value),
        _reference(reference)
    {}

    ParameterComponent(ParameterComponent const & other) :
        Base(other._value),
        _reference(other._reference)
    {}

    lsst::afw::geom::Point2D const getPosition() const { return _reference + _value; }

    lsst::afw::geom::Point2D const & getReference() const { return _reference; }
    
protected:

    template <ParameterType E> friend class grid::ParameterComponent;

    void _writeParameters(double * paramIter) const {
        paramIter[0] = _value.getX();
        paramIter[1] = _value.getY();
    }

    void _readParameters(double const * paramIter) {
        _value.setX(paramIter[0]);
        _value.setY(paramIter[1]); 
    }

    lsst::afw::geom::Point2D _reference;
};

typedef ParameterComponent<POSITION> PositionComponent;

/**
 *  @brief A radius parameter.
 */
template<>
class ParameterComponent<RADIUS> 
    : public ParameterComponentBase<ParameterComponent<RADIUS>,Radius,1>
{
    typedef ParameterComponentBase<ParameterComponent<RADIUS>,Radius,1> Base;
public:

    explicit ParameterComponent(Value const & value) : Base(value) {}

    ParameterComponent(ParameterComponent const & other) : Base(other._value) {}

protected:

    template <ParameterType E> friend class grid::ParameterComponent;

    void _writeParameters(double * paramIter) const { paramIter[0] = _value; }

    void _readParameters(double const * paramIter) { _value = paramIter[0]; }

};

typedef ParameterComponent<RADIUS> RadiusComponent;


/**
 *  @brief A pair of ellipticity parameters.
 */
template <>
class ParameterComponent<ELLIPTICITY>
    : public ParameterComponentBase<ParameterComponent<ELLIPTICITY>,Ellipticity,2> 
{
    typedef ParameterComponentBase<ParameterComponent<ELLIPTICITY>, Ellipticity,2> Base;
public:

    explicit ParameterComponent(Value const & value) : Base(value) {}

    ParameterComponent(ParameterComponent const & other) : Base(other._value) {}

protected:

    template <ParameterType E> friend class grid::ParameterComponent;

    void _writeParameters(double * paramIter) const {
        paramIter[0] = _value.getE1();
        paramIter[1] = _value.getE2();
    }

    void _readParameters(double const * paramIter) {
        _value.setE1(paramIter[0]);
        _value.setE2(paramIter[1]);
    }
};

typedef ParameterComponent<ELLIPTICITY> EllipticityComponent;

}}}} // namespace lsst::meas::multifit::definition

#endif // !LSST_MEAS_MULTIFIT_DEFINITION_parameters
