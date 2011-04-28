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

namespace detail {

template <ParameterType E> struct ParameterComponentTraits;

template <>
struct ParameterComponentTraits<POSITION> {

    typedef lsst::afw::geom::Point2D Value;
    static int const SIZE = 2;

    static void readParameters(double const * paramIter, Value & value) {
        value[0] += paramIter[0];
        value[1] += paramIter[1];
    }
    static void writeParameters(double * paramIter, Value const & value) {
        paramIter[0] = 0.0;
        paramIter[1] = 0.0;
    }
    static void printValue(std::ostream & os, Value const & value) { os << value; }
};

template <>
struct ParameterComponentTraits<RADIUS> {

    typedef Radius Value;
    static int const SIZE = 1;

    static void readParameters(double const * paramIter, Value & value) {
        value = paramIter[0];
    }
    static void writeParameters(double * paramIter, Value const & value) {
        paramIter[0] = value;
    }
    static void printValue(std::ostream & os, Value const & value) { os << (double)value; }

};

template <>
struct ParameterComponentTraits<ELLIPTICITY> {

    typedef Ellipticity Value;
    static int const SIZE = 2;

    static void readParameters(double const * paramIter, Value & value) {
        value.setE1(paramIter[0]);
        value.setE2(paramIter[1]);
    }
    static void writeParameters(double * paramIter, Value const & value) {
        paramIter[0] = value.getE1();
        paramIter[1] = value.getE2();
    }
    static void printValue(std::ostream & os, Value const & value) { os << value.getComplex(); }

};

template <ParameterType E>
class ParameterComponentBase {
public:

    typedef typename detail::ParameterComponentTraits<E>::Value Value;
    static int const SIZE = detail::ParameterComponentTraits<E>::SIZE;

    /**
     *  @brief True if the position is allowed to vary during fitting.
     */
    bool const isActive() const { return _active; }

    /**
     *  @brief Return the value of the ParameterComponent.
     *
     *  For most parameters, this is the same as the values in the initial parameter vector.
     *  For POSITION, the parameter vector contains a pair of (x,y) offsets from this value.
     */
    Value const & getValue() const { return _value; }

protected:

    explicit ParameterComponentBase(Value const & value, bool active) :
        _active(active), _value(value)
    {}

    ParameterComponentBase(ParameterComponentBase const & other) :
        _active(other._active), _value(other._value)
    {}

    bool _active;
    Value _value;
};

} // namespace detail

namespace definition {

template <ParameterType E>
class ParameterComponent : public detail::ParameterComponentBase<E> {
public:
    
    typedef boost::shared_ptr< ParameterComponent<E> > Ptr;
    typedef boost::shared_ptr< ParameterComponent<E> const > ConstPtr;
    typedef typename detail::ParameterComponentTraits<E>::Value Value;

#ifndef SWIG // these are wrapped explicitly; SWIG is confused by the typedefs and "bool &"

    //@{
    /**
     *  @brief True if the position is allowed to vary during fitting.
     */
    bool & isActive() { return this->_active; }
    void setActive(bool active) { this->_active = active; }
    //@}
    
    //@{
    /**
     *  @brief Return and/or set the value of the ParameterComponent.
     *
     *  For most parameters, this is the same as the values in the initial parameter vector.
     *  For POSITION, the parameter vector contains a pair of (x,y) offsets from this value.
     *
     *  In Python, this is additionally wrapped as a property 'component.value',
     *  because the syntax 'component.getValue() = ' is illegal.
     */
    Value & getValue() { return this->_value; }
    void setValue(Value const & value) { this->_value = value; }
    //@}

    // Use const accessors from base class.
    using detail::ParameterComponentBase<E>::getValue;
    using detail::ParameterComponentBase<E>::isActive;

    /**
     *  @brief Deep-copy.
     *
     *  @note Not called clone() because it's not virtual, and it doesn't need to be.
     */
    Ptr copy() const { return Ptr(new ParameterComponent(*this)); }

    /**
     *  @brief Create a new ParameterComponent.
     *
     *  Constructors are private to ensure we only get shared_ptrs to these things.
     */
    static Ptr make(Value const & value, bool active=true) {
        return Ptr(new ParameterComponent(value, active));
    }

#endif

private:

    ParameterComponent(ParameterComponent const & other) : detail::ParameterComponentBase<E>(*this) {}

    explicit ParameterComponent(Value const & value, bool active) : 
        detail::ParameterComponentBase<E>(value, active) {}

};

typedef ParameterComponent<POSITION> PositionComponent;
typedef ParameterComponent<RADIUS> RadiusComponent;
typedef ParameterComponent<ELLIPTICITY> EllipticityComponent;

template <ParameterType E>
std::ostream & operator<<(std::ostream & os, ParameterComponent<E> const & component);

}}}} // namespace lsst::meas::multifit::definition

#endif // !LSST_MEAS_MULTIFIT_DEFINITION_parameters
