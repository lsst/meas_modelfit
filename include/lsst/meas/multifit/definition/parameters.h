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

template <ParameterType P> class ParameterComponent;

/**
 *  @brief Parameters representing a position.
 *
 *  A position is represented as an offset from a reference point, which is stored within the
 *  position.
 */
template <>
class ParameterComponent<POSITION> {
public:
    
    typedef boost::shared_ptr< ParameterComponent<POSITION> > Ptr;
    typedef boost::shared_ptr< ParameterComponent<POSITION> const > ConstPtr;

    typedef lsst::afw::geom::Extent2D Value;
    static int const SIZE = 2;

    explicit ParameterComponent(lsst::afw::geom::Point2D const & reference, Value const & value = Value()) :
        active(true), _value(value), _reference(reference)
    {}

    ParameterComponent(ParameterComponent const & other) :
        active(other.active), _value(other._value), _reference(other._reference)
    {}

    /// @brief True if the position is allowed to vary during fitting.
    bool active;

    /**
     *  @brief Initial value of the parameter.
     *
     *  For position, this is generally zero because the parameters are the offsets from the reference
     *  position.
     */
    Value const & getValue() const { return _value; }

    /// @brief Deep-copy.  Not clone() because it's not virtual, and doesn't need to be.
    Ptr copy() const { return boost::make_shared< ParameterComponent<POSITION> >(*this); }

    /// @brief Return the initial value of the position.
    lsst::afw::geom::Point2D const getPosition() const { return _reference + _value; }

    /// @brief Return the reference posiion.x
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

    Value _value;
    lsst::afw::geom::Point2D _reference;
};

typedef ParameterComponent<POSITION> PositionComponent;

std::ostream & operator<<(std::ostream & os, PositionComponent const & component);

/**
 *  @brief A radius parameter.
 */
template<>
class ParameterComponent<RADIUS> {
public:

    typedef boost::shared_ptr< ParameterComponent<RADIUS> > Ptr;
    typedef boost::shared_ptr< ParameterComponent<RADIUS> const > ConstPtr;

    typedef Radius Value;
    static int const SIZE = 1;

    explicit ParameterComponent(Value const & value) : active(true), _value(value) {}

    ParameterComponent(ParameterComponent const & other) : active(other.active), _value(other._value) {}

    /// @brief True if the radius is allowed to vary during fitting.
    bool active;

    /// @brief Initial value of the parameter.
    Value const & getValue() const { return _value; }

    /// @brief Deep-copy.  Not clone() because it's not virtual, and doesn't need to be.
    Ptr copy() const { return boost::make_shared< ParameterComponent<RADIUS> >(*this); }

protected:

    template <ParameterType E> friend class grid::ParameterComponent;

    void _writeParameters(double * paramIter) const { paramIter[0] = _value; }

    void _readParameters(double const * paramIter) { _value = paramIter[0]; }

    Value _value;
};

typedef ParameterComponent<RADIUS> RadiusComponent;

std::ostream & operator<<(std::ostream & os, RadiusComponent const & component);

/**
 *  @brief A pair of ellipticity parameters.
 */
template <>
class ParameterComponent<ELLIPTICITY> {
public:

    typedef boost::shared_ptr< ParameterComponent<ELLIPTICITY> > Ptr;
    typedef boost::shared_ptr< ParameterComponent<ELLIPTICITY> const > ConstPtr;

    typedef Ellipticity Value;
    static int const SIZE = 2;

    explicit ParameterComponent(Value const & value) : active(true), _value(value) {}

    ParameterComponent(ParameterComponent const & other) : active(other.active), _value(other._value) {}

    /// @brief True if the ellipticity is allowed to vary during fitting.
    bool active;

    /// @brief Initial value of the parameter.
    Value const & getValue() const { return _value; }

    /// @brief Deep-copy.  Not clone() because it's not virtual, and doesn't need to be.
    Ptr copy() const { return boost::make_shared< ParameterComponent<ELLIPTICITY> >(*this); }

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

    Value _value;
};

typedef ParameterComponent<ELLIPTICITY> EllipticityComponent;

std::ostream & operator<<(std::ostream & os, EllipticityComponent const & component);

}}}} // namespace lsst::meas::multifit::definition

#endif // !LSST_MEAS_MULTIFIT_DEFINITION_parameters
