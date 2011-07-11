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

#ifndef LSST_MEAS_MULTIFIT_DEFINITION_SharedElement
#define LSST_MEAS_MULTIFIT_DEFINITION_SharedElement

#include "lsst/meas/multifit/constants.h"
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

namespace lsst { namespace meas { namespace multifit {

namespace detail {

template <SharedElementType E> struct SharedElementTraits;

struct CircleConstraint {
    double max;

    explicit CircleConstraint(double max_) : max(max_) {}

    friend inline std::ostream & operator<<(std::ostream & os, CircleConstraint const & k) {
        return os << "||.|| <  " << k.max;
    }

private:
    template <SharedElementType E> friend class grid::SharedElement;

    /// Return true if the parameters are in-bounds.
    bool checkBounds(double const * parameterIter) const {
        return parameterIter[0] * parameterIter[0] + parameterIter[1] + parameterIter[1] <= max * max;
    }

    /**
     *  If the parameters are out of bounds, move them to the boundary and return
     *  a positive value that increases as the necessary parameter change increases.
     *  Return 0.0 if the parameters are already in-bounds
     */
    double clipToBounds(double * parameterIter) const {
        double a = std::sqrt(max * max / (parameterIter[0] * parameterIter[0] + parameterIter[1] + parameterIter[1]));
        if (a < 1.0) {
            parameterIter[0] *= a;
            parameterIter[1] *= a;
            return -std::log(a);
        }
        return 0.0;
    }
};

struct MinMaxConstraint {
    double min;
    double max;

    /// Return true if the parameters are in-bounds.
    explicit MinMaxConstraint(double min_=0.0, double max_=std::numeric_limits<double>::infinity()) 
        : min(min_), max(max_) {}

    friend inline std::ostream & operator<<(std::ostream & os, MinMaxConstraint const & k) {
        return os << k.min << " <= (.) <= " << k.max;
    }

private:
    template <SharedElementType E> friend class grid::SharedElement;

    /// Return true if the parameters are in-bounds.
    bool checkBounds(double const * parameterIter) const {
        return *parameterIter >= min && *parameterIter <= max;
    }

    /**
     *  If the parameters are out of bounds, move them to the boundary and return
     *  a positive value that increases as the necessary parameter change increases.
     *  Return 0.0 if the parameters are already in-bounds
     */
    double clipToBounds(double * parameterIter) const {
        if (*parameterIter < min) {
            double r = min - *parameterIter;
            *parameterIter = min;
            return r;
        } else if (*parameterIter > max) {
            double r = *parameterIter - max;
            *parameterIter = max;
            return r;
        }
        return 0.0;
    }
};

template <>
struct SharedElementTraits<POSITION> {

    typedef lsst::afw::geom::Point2D Value;
    typedef CircleConstraint Bounds;
    static int const SIZE = 2;

    static void readParameters(double const * parameterIter, Value & value) {
        value[0] += parameterIter[0];
        value[1] += parameterIter[1];
    }
    static void writeParameters(double * parameterIter, Value const & value) {
        parameterIter[0] = 0.0;
        parameterIter[1] = 0.0;
    }
    static void printValue(std::ostream & os, Value const & value) { os << value; }

    static Bounds getDefaultBounds() { return Bounds(5.0); }
};

template <>
struct SharedElementTraits<RADIUS> {

    typedef Radius Value;
    typedef MinMaxConstraint Bounds;
    static int const SIZE = 1;

    static void readParameters(double const * parameterIter, Value & value) {
        value = parameterIter[0];
    }
    static void writeParameters(double * parameterIter, Value const & value) {
        parameterIter[0] = value;
    }
    static void printValue(std::ostream & os, Value const & value) { os << (double)value; }

    static Bounds getDefaultBounds() {
        return Bounds(std::numeric_limits<double>::epsilon(), std::numeric_limits<double>::infinity());
    }
};

template <>
struct SharedElementTraits<ELLIPTICITY> {

    typedef Ellipticity Value;
    typedef CircleConstraint Bounds;
    static int const SIZE = 2;

    static void readParameters(double const * parameterIter, Value & value) {
        value.setE1(parameterIter[0]);
        value.setE2(parameterIter[1]);
    }
    static void writeParameters(double * parameterIter, Value const & value) {
        parameterIter[0] = value.getE1();
        parameterIter[1] = value.getE2();
    }
    static void printValue(std::ostream & os, Value const & value) { os << value.getComplex(); }

    static Bounds getDefaultBounds() {
        return Bounds(-std::log(std::numeric_limits<double>::epsilon()));
    }
};

template <SharedElementType E>
class SharedElementBase {
public:

    typedef typename detail::SharedElementTraits<E>::Value Value;
    typedef typename detail::SharedElementTraits<E>::Bounds Bounds;
    static int const SIZE = detail::SharedElementTraits<E>::SIZE;

    /**
     *  @brief True if the position is allowed to vary during fitting.
     */
    bool const isActive() const { return _active; }

    /**
     *  @brief Return the value of the SharedElement.
     *
     *  For most parameters, this is the same as the values in the initial parameter vector.
     *  For POSITION, the parameter vector contains a pair of (x,y) offsets from this value.
     */
    Value const & getValue() const { return _value; }

    Bounds const & getBounds() const { return _bounds; }

protected:

    explicit SharedElementBase(Value const & value, Bounds const & bounds, bool active) :
        _active(active), _value(value), _bounds(bounds)
    {}

    SharedElementBase(SharedElementBase const & other) :
        _active(other._active), _value(other._value), _bounds(other._bounds)
    {}

    bool _active;
    Value _value;
    Bounds _bounds;
};

} // namespace detail

namespace definition {

template <SharedElementType E>
class SharedElement : public detail::SharedElementBase<E> {
public:
    
    typedef boost::shared_ptr< SharedElement<E> > Ptr;
    typedef boost::shared_ptr< SharedElement<E> const > ConstPtr;
    typedef typename detail::SharedElementTraits<E>::Value Value;
    typedef typename detail::SharedElementTraits<E>::Bounds Bounds;

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
     *  @brief Return and/or set the value of the SharedElement.
     *
     *  For most parameters, this is the same as the values in the initial parameter vector.
     *  For POSITION, the parameter vector contains a pair of (x,y) offsets from this value.
     */
    Value & getValue() { return this->_value; }
    void setValue(Value const & value) { this->_value = value; }
    //@}
    
    //@{
    /**
     *  @brief Return and/or set the bounds of the SharedElement.
     */
    Bounds & getBounds() { return this->_bounds; }
    void setBounds(Bounds const & bounds) { this->_bounds = bounds; }
    //@}

    // Use const accessors from base class.
    using detail::SharedElementBase<E>::getValue;
    using detail::SharedElementBase<E>::getBounds;
    using detail::SharedElementBase<E>::isActive;

    /**
     *  @brief Deep-copy.
     *
     *  @note Not called clone() because it's not virtual, and it doesn't need to be.
     */
    Ptr copy() const { return Ptr(new SharedElement(*this)); }

    /**
     *  @brief Create a new SharedElement.
     *
     *  Constructors are private to ensure we only get shared_ptrs to these things.
     */
    static Ptr make(Value const & value, bool active=true) {
        return Ptr(
            new SharedElement(
                value, detail::SharedElementTraits<E>::getDefaultBounds(), active
            )
        );
    }

    /**
     *  @brief Create a new SharedElement.
     *
     *  Constructors are private to ensure we only get shared_ptrs to these things.
     */
    static Ptr make(Value const & value, Bounds const & bounds, bool active=true) {
        return Ptr(new SharedElement(value, bounds, active));
    }

#endif

private:

    SharedElement(SharedElement const & other) : detail::SharedElementBase<E>(other) {}

    explicit SharedElement(Value const & value, Bounds const & bounds, bool active) : 
        detail::SharedElementBase<E>(value, bounds, active) {}

};

typedef SharedElement<POSITION> PositionElement;
typedef SharedElement<RADIUS> RadiusElement;
typedef SharedElement<ELLIPTICITY> EllipticityElement;

#ifndef SWIG
template <SharedElementType E>
std::ostream & operator<<(std::ostream & os, SharedElement<E> const & component);
#endif

}}}} // namespace lsst::meas::multifit::definition

#endif // !LSST_MEAS_MULTIFIT_DEFINITION_SharedElement
