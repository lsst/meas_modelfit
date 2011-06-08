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

#ifndef LSST_MEAS_MULTIFIT_DEFINITION_ObjectComponent
#define LSST_MEAS_MULTIFIT_DEFINITION_ObjectComponent


#include "lsst/meas/multifit/definition/elements.h"
#include "lsst/meas/multifit/ModelBasis.h"

namespace lsst { namespace meas { namespace multifit {

namespace detail {

/// @brief Utility (nonpolymorphic) base class for definition::ObjectComponent and grid::ObjectComponent.
class ObjectComponentBase {
public:

    ID const id;

    /// @brief Return the basis used for extended objects.
    ModelBasis::Ptr const getBasis() const { return _basis; }

    /// @brief Return the ratio of the actual radius of this object to the radius parameter.
    double const getRadiusFactor() const { return _radiusFactor; }

    /**
     *  @brief Return whether the object is time-variable (has different coefficients
     *         in each exposure).
     */
    bool const isVariable() const { return _variable; };

protected:
    
    ObjectComponentBase(
        ID const id_, 
        double const radiusFactor,
        bool const variable,
        ModelBasis::Ptr const & basis = ModelBasis::Ptr()
    ) : id(id_), _radiusFactor(radiusFactor), _variable(variable), _basis(basis) {}

    ObjectComponentBase(ObjectComponentBase const & other) : 
        id(other.id), _radiusFactor(other._radiusFactor), _variable(other._variable), _basis(other._basis)
    {}

    double _radiusFactor;
    bool _variable;
    ModelBasis::Ptr _basis;
};

} // namespace detail


namespace definition {

/**
 *  @brief An astronomical object in a multifit definition.
 *
 *  An ObjectComponent will typically represent a point source or elliptical source.  For the former, the
 *  radius, ellipticity, and basis pointers will be empty, while they will all be non-empty for the latter.
 *
 *  Most of the accessors of object return by non-const reference, reflecting that ObjectComponent is essentially
 *  a struct that is validated only when a Grid is constructed from the Definition.  We have used accessors
 *  rather than public data members so the interface is similar to that of grid::ObjectComponent, which behaves
 *  like a const version of definition::ObjectComponent.
 *
 *  In Python, accessors of ObjectComponent will be wrapped both with "get/set" methods and as Python properties.
 */
class ObjectComponent : public detail::ObjectComponentBase {
public:

    explicit ObjectComponent(ID id_) : detail::ObjectComponentBase(id_, 1.0, false) {}

    ObjectComponent(ObjectComponent const & other) :
        detail::ObjectComponentBase(other),
        _position(other._position), _radius(other._radius), _ellipticity(other._ellipticity)
    {}

    static ObjectComponent makeStar(
        ID const id, 
        lsst::afw::geom::Point2D const & position, 
        bool const isVariable = false,
        bool const isPositionActive=false
    );

    static ObjectComponent makeGalaxy(
        ID const id,
        ModelBasis::Ptr const & basis,
        lsst::afw::geom::ellipses::Ellipse const & ellipse,
        bool const isEllipticityActive=false,
        bool const isRadiusActive=false,
        bool const isPositionActive=false
    );

    //@{
    /**
     *  @brief Return and/or set the parameter components.
     */
#ifndef SWIG
    PositionElement::Ptr & getPosition() { return _position; }
    RadiusElement::Ptr & getRadius() { return _radius; }
    EllipticityElement::Ptr & getEllipticity() { return _ellipticity; }

    PositionElement::Ptr const getPosition() const { return _position; }
    RadiusElement::Ptr const getRadius() const { return _radius; }
    EllipticityElement::Ptr const getEllipticity() const { return _ellipticity; }

    template <ParameterType E>
    typename ParameterElement<E>::Ptr & getElement() {
        typename ParameterElement<E>::Ptr const * p;
        getElementImpl(p);
        return const_cast<typename ParameterElement<E>::Ptr &>(*p);
    }

    template <ParameterType E>
    typename ParameterElement<E>::Ptr const getElement() const {
        typename ParameterElement<E>::Ptr const * p;
        getElementImpl(p);
        return *p;
    }
    //@}

    /// @brief Return and/or set the basis used for extended objects.
    ModelBasis::Ptr & getBasis() { return _basis; }

    /// @brief Return and/or set the ratio of the actual radius of this object to the radius parameter.
    double & getRadiusFactor() { return _radiusFactor; }

    /**
     *  @brief Return and/or set whether the object is time-variable (has different coefficients
     *         in each exposure).
     */
    bool & isVariable() { return _variable; };

#endif

    using detail::ObjectComponentBase::getBasis;
    using detail::ObjectComponentBase::getRadiusFactor;
    using detail::ObjectComponentBase::isVariable;

    //@{
    /**
     *  @brief Setters for basis, radius factor and variability.
     *
     *  These aren't necessary in C++ because the non-const getters return references, but they might
     *  be expected, and they're the only way to do things from Python.
     */
    void setBasis(ModelBasis::Ptr const & basis) { _basis = basis; }
    void setRadiusFactor(double factor) { _radiusFactor = factor; }
    void setVariable(bool variable) { _variable = variable; }
    //@}

private:

    friend class grid::Initializer;

    explicit ObjectComponent(detail::ObjectComponentBase const & other) : detail::ObjectComponentBase(other) {}

    void getElementImpl(PositionElement::Ptr const * & p) const { p = &_position; }
    void getElementImpl(RadiusElement::Ptr const * & p) const { p = &_radius; }
    void getElementImpl(EllipticityElement::Ptr const * & p) const { p = &_ellipticity; }

    PositionElement::Ptr _position;
    RadiusElement::Ptr _radius;
    EllipticityElement::Ptr _ellipticity;
};

#ifndef SWIG
std::ostream & operator<<(std::ostream & os, ObjectComponent const & obj);
#endif

}}}} // namespace lsst::meas::multifit::definition

#endif // !LSST_MEAS_MULTIFIT_DEFINITION_ObjectComponent
