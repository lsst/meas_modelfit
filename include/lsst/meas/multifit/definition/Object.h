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

#ifndef LSST_MEAS_MULTIFIT_DEFINITION_Object
#define LSST_MEAS_MULTIFIT_DEFINITION_Object

#include "multifit/definition/parameters.hpp"
#include <modeling/EllipseBasis.hpp>
#include <modeling/LinearConstraint.hpp>

namespace lsst { namespace meas { namespace multifit { namespace definition {

class InvalidDefinitionError : public std::invalid_argument {
public:
    InvalidDefinitionError(char const * msg) : std::invalid_argument(msg) {}
};

class Object {
public:

    explicit Object(ID id_) : id(id_), radius_factor(1.0), is_variable(false) {}

    Object(Object const & other) : 
        id(other.id), position(other.position),
        radius(other.radius), ellipticity(other.ellipticity),
        basis(other.basis),
        radius_factor(other.radius_factor),
        is_variable(other.is_variable) 
    {}

    static Object makeStar(
        ID id, 
        agl::PointD const & position, 
        bool is_variable = false
    );

    static Object makeGalaxy(
        ID id,
        modeling::EllipseBasis::Ptr const & basis,
        agl::Ellipse const & ellipse
    );        

    ID const id;

    PositionComponent::Ptr position;
    RadiusComponent::Ptr radius;
    EllipticityComponent::Ptr ellipticity;

    modeling::EllipseBasis::Ptr basis;
    double radius_factor;
    bool is_variable;

    modeling::LinearConstraint::Ptr inequality_constraint;
    modeling::LinearConstraint::Ptr equality_constraint;

};

}}}} // namespace lsst::meas::multifit::definition

#endif // !LSST_MEAS_MULTIFIT_DEFINITION_Object
