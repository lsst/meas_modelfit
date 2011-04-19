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

namespace lsst { namespace meas { namespace multifit { namespace definition {

LSST_EXCEPTION_TYPE(InvalidDefinitionError,
                    lsst::pex::exceptions::InvalidParameterException,
                    lsst::meas::multifit::definition::InvalidDefinitionError);

class Object {
public:

    explicit Object(ID id_) : id(id_), radiusFactor(1.0), isVariable(false) {}

    Object(Object const & other) : 
        id(other.id), position(other.position),
        radius(other.radius), ellipticity(other.ellipticity),
        basis(other.basis),
        radiusFactor(other.radiusFactor),
        isVariable(other.isVariable) 
    {}

    static Object makeStar(
        ID id, 
        lsst::afw::geom::Point2D const & position, 
        bool isVariable = false,
        bool isPositionActive=false
    );

    static Object makeGalaxy(
        ID id,
        ModelBasis::Ptr const & basis,
        lsst::afw::geom::ellipses::Ellipse const & ellipse,
        bool isEllipticityActive=false,
        bool isRadiusActive=false,
        bool isPositionActive=false
    );        

    /// @brief Throw an exception (LogicErrorException) if the object lacks a radius or ellipticity.
    void requireEllipse() const;

    ID const id;

    PositionComponent::Ptr position;
    RadiusComponent::Ptr radius;
    EllipticityComponent::Ptr ellipticity;

    ModelBasis::Ptr basis;
    double radiusFactor;
    bool isVariable;

};

}}}} // namespace lsst::meas::multifit::definition

#endif // !LSST_MEAS_MULTIFIT_DEFINITION_Object
