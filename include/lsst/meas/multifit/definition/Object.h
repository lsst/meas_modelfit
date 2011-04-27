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

#ifndef LSST_MEAS_MULTIFIT_DEFINITION_Object
#define LSST_MEAS_MULTIFIT_DEFINITION_Object


#include "lsst/meas/multifit/definition/parameters.h"
#include "lsst/meas/multifit/ModelBasis.h"

namespace lsst { namespace meas { namespace multifit {

namespace detail {

/// @brief Utility (nonpolymorphic) base class for definition::Object and grid::Object.
class ObjectBase {
public:

    ID const id;

    /// @brief Return the basis used for extended objects.
    ModelBasis::ConstPtr const getBasis() const { return _basis; }

    /// @brief Return the ratio of the actual radius of this object to the radius parameter.
    double const getRadiusFactor() const { return _radiusFactor; }

    /**
     *  @brief Return whether the object is time-variable (has different coefficients
     *         in each exposure).
     */
    bool const isVariable() const { return _variable; };

protected:
    
    ObjectBase(
        ID id_, 
        double radiusFactor,
        bool variable,
        ModelBasis::ConstPtr const & basis = ModelBasis::ConstPtr()
    ) : id(id_), _radiusFactor(radiusFactor), _variable(variable), _basis(basis) {}

    ObjectBase(ObjectBase const & other) : 
        id(other.id), _radiusFactor(other._radiusFactor), _variable(other._variable), _basis(other._basis)
    {}

    double _radiusFactor;
    bool _variable;
    ModelBasis::ConstPtr _basis;
};

} // namespace detail


namespace definition {

/**
 *  @brief An astronomical object in a multifit definition.
 *
 *  An Object will typically represent a point source or elliptical source.  For the former, the
 *  radius, ellipticity, and basis pointers will be empty, while they will all be non-empty for the latter.
 *
 *  Most of the accessors of object return by non-const reference, reflecting that Object is essentially
 *  a struct that is validated only when a Grid is constructed from the Definition.  We have used accessors
 *  rather than public data members so the interface is similar to that of grid::Object, which behaves
 *  like a const version of definition::Object.
 *
 *  In Python, accessors of Object will be wrapped both with "get/set" methods and as Python properties.
 */
class Object : public detail::ObjectBase {
public:

    explicit Object(ID id_) : detail::ObjectBase(id_, 1.0, false) {}

    Object(Object const & other) :
        detail::ObjectBase(other),
        _position(other._position), _radius(other._radius), _ellipticity(other._ellipticity)
    {}

    static Object makeStar(
        ID id, 
        lsst::afw::geom::Point2D const & position, 
        bool isVariable = false,
        bool isPositionActive=false
    );

    static Object makeGalaxy(
        ID id,
        ModelBasis::ConstPtr const & basis,
        lsst::afw::geom::ellipses::Ellipse const & ellipse,
        bool isEllipticityActive=false,
        bool isRadiusActive=false,
        bool isPositionActive=false
    );

    //@{
    /**
     *  @brief Return and/or set the parameter components.
     */
    PositionComponent::Ptr & getPosition() { return _position; }
    RadiusComponent::Ptr & getRadius() { return _radius; }
    EllipticityComponent::Ptr & getEllipticity() { return _ellipticity; }

    PositionComponent::Ptr const getPosition() const { return _position; }
    RadiusComponent::Ptr const getRadius() const { return _radius; }
    EllipticityComponent::Ptr const getEllipticity() const { return _ellipticity; }

    template <ParameterType E>
    typename ParameterComponent<E>::Ptr & getComponent() {
        typename ParameterComponent<E>::Ptr const * p;
        getComponentImpl(p);
        return const_cast<typename ParameterComponent<E>::Ptr &>(*p);
    }

    template <ParameterType E>
    typename ParameterComponent<E>::Ptr const getComponent() const {
        typename ParameterComponent<E>::Ptr const * p;
        getComponentImpl(p);
        return *p;
    }
    //@}

    /// @brief Return and/or set the basis used for extended objects.
    ModelBasis::ConstPtr & getBasis() { return _basis; }

    /// @brief Return and/or set the ratio of the actual radius of this object to the radius parameter.
    double & getRadiusFactor() { return _radiusFactor; }

    /**
     *  @brief Return and/or set whether the object is time-variable (has different coefficients
     *         in each exposure).
     */
    bool & isVariable() { return _variable; };

    using detail::ObjectBase::getBasis;
    using detail::ObjectBase::getRadiusFactor;
    using detail::ObjectBase::isVariable;

private:

    friend class grid::Initializer;

    explicit Object(detail::ObjectBase const & other) : detail::ObjectBase(other) {}

    void getComponentImpl(PositionComponent::Ptr const * & p) const { p = &_position; }
    void getComponentImpl(RadiusComponent::Ptr const * & p) const { p = &_radius; }
    void getComponentImpl(EllipticityComponent::Ptr const * & p) const { p = &_ellipticity; }

    PositionComponent::Ptr _position;
    RadiusComponent::Ptr _radius;
    EllipticityComponent::Ptr _ellipticity;
};

std::ostream & operator<<(std::ostream & os, Object const & obj);

}}}} // namespace lsst::meas::multifit::definition

#endif // !LSST_MEAS_MULTIFIT_DEFINITION_Object
