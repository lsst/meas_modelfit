// -*- lsst-c++ -*-

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

%define %PointerEQ(CLS)
%extend CLS {
    %feature("shadow") _equals %{
        def __eq__(self, other):
            try:
                return $action(self, other)
            except ArgumentError:
                return NotImplemented
    %}
    %pythoncode %{
        def __ne__(self, other):
            return not self.__eq__(other)
    %}
     bool _equals(CLS * other) {
        return self == other;
    }
}
%enddef

%define %AddStreamRepr(CLS)
%extend CLS {
    std::string __repr__() {
        std::ostringstream os;
        os << (*self);
        return os.str();
    }
}
%enddef

//----------------------------- SharedElements -----------------------------------

%{
#include "lsst/meas/multifit/definition/SharedElement.h"
%}

SWIG_SHARED_PTR(definition_PositionElementPtr, lsst::meas::multifit::definition::SharedElement<lsst::meas::multifit::POSITION>);
SWIG_SHARED_PTR(definition_RadiusElementPtr, lsst::meas::multifit::definition::SharedElement<lsst::meas::multifit::RADIUS>);
SWIG_SHARED_PTR(definition_EllipticityElementPtr, lsst::meas::multifit::definition::SharedElement<lsst::meas::multifit::ELLIPTICITY>);

%define %DeclareDefinitionSharedElement(TITLE, LOWER, UPPER, CONSTRAINT)
%template(definition_##TITLE##Element)
lsst::meas::multifit::definition::SharedElement<lsst::meas::multifit::UPPER>;
%PointerEQ(lsst::meas::multifit::definition::SharedElement<lsst::meas::multifit::UPPER>)
%AddStreamRepr(lsst::meas::multifit::definition::SharedElement<lsst::meas::multifit::UPPER>)
%extend lsst::meas::multifit::definition::SharedElement<lsst::meas::multifit::UPPER> {
    static boost::shared_ptr< 
        lsst::meas::multifit::definition::SharedElement< lsst::meas::multifit::UPPER >
    > make(lsst::meas::multifit::TITLE const & value, bool active=true) {
        return lsst::meas::multifit::definition::SharedElement< lsst::meas::multifit::UPPER >::make(
            value, active
        );
    }
    lsst::meas::multifit::TITLE & getValue() {
        return self->getValue();
    }
    void setValue(lsst::meas::multifit::TITLE const & value) {
        self->getValue() = value;
    }
    lsst::meas::multifit::CONSTRAINT & getBounds() {
        return self->getBounds();
    }
    void setBounds(lsst::meas::multifit::CONSTRAINT const & bounds) {
        self->getBounds() = bounds;
    }
    bool isActive() {
        return self->isActive();
    }
    void setActive(bool active) {
        self->isActive() = active;
    }
}
%extend lsst::meas::multifit::definition::ObjectComponent {
    boost::shared_ptr< lsst::meas::multifit::definition::SharedElement< lsst::meas::multifit::UPPER > > get ## TITLE() {
        return self->get##TITLE();
    }
    void set ## TITLE(boost::shared_ptr< lsst::meas::multifit::definition::SharedElement< lsst::meas::multifit::UPPER > > value) {
        self->get ## TITLE() = value;
    }
}
%enddef

%rename(detail_CircleConstraint) lsst::meas::multifit::detail::CircleConstraint;
%rename(detail_MinMaxConstraint) lsst::meas::multifit::detail::MinMaxConstraint;

%include "lsst/meas/multifit/definition/SharedElement.h"

//------------------------------------- ObjectComponent ---------------------------------------

%{
#include "lsst/meas/multifit/definition/ObjectComponent.h"
%}

SWIG_SHARED_PTR(detail_ObjectComponentBasePtr, lsst::meas::multifit::detail::ObjectComponentBase);
SWIG_SHARED_PTR_DERIVED(definition_ObjectComponentPtr, lsst::meas::multifit::detail::ObjectComponentBase, lsst::meas::multifit::definition::ObjectComponent);

%rename(definition_ObjectComponent) lsst::meas::multifit::definition::ObjectComponent;

%include "lsst/meas/multifit/definition/ObjectComponent.h"

%DeclareDefinitionSharedElement(Position, position, POSITION, detail::CircleConstraint);
%DeclareDefinitionSharedElement(Radius, radius, RADIUS, detail::MinMaxConstraint);
%DeclareDefinitionSharedElement(Ellipticity, ellipticity, ELLIPTICITY, detail::CircleConstraint);

%PointerEQ(lsst::meas::multifit::definition::ObjectComponent)
%AddStreamRepr(lsst::meas::multifit::definition::ObjectComponent)

//------------------------------------- Frame ---------------------------------------

%{
#include "lsst/meas/multifit/definition/Frame.h"
%}

SWIG_SHARED_PTR(detail_FrameBasePtr, lsst::meas::multifit::detail::FrameBase);
SWIG_SHARED_PTR_DERIVED(definition_FramePtr, lsst::meas::multifit::detail::FrameBase, lsst::meas::multifit::definition::Frame);

%rename(definition_Frame) lsst::meas::multifit::definition::Frame;

%include "lsst/meas/multifit/definition/Frame.h"

%PointerEQ(lsst::meas::multifit::definition::Frame)
%AddStreamRepr(lsst::meas::multifit::definition::Frame)

%extend lsst::meas::multifit::definition::Frame {
    lsst::meas::multifit::FilterId getFilterId() { return self->getFilterId(); }
    lsst::meas::multifit::Wcs::Ptr getWcs() { return self->getWcs(); }
    lsst::meas::multifit::Psf::Ptr getPsf() { return self->getPsf(); }
    lsst::meas::multifit::Footprint::Ptr getFootprint() { return self->getFootprint(); }
    lsst::ndarray::Array<lsst::meas::multifit::Pixel,1,1> getData() { return self->getData(); }
    lsst::ndarray::Array<lsst::meas::multifit::Pixel,1,1> getWeights() { return self->getWeights(); }
}

%template(make) lsst::meas::multifit::definition::Frame::make<float>;
%template(make) lsst::meas::multifit::definition::Frame::make<double>;


//-------------------------------------- Set -----------------------------------------

%include "lsst/meas/multifit/containers.i"

%DeclareMutableSet(lsst::meas::multifit::definition::ObjectComponent, definition_ObjectComponentSet)
%DeclareMutableSet(lsst::meas::multifit::definition::Frame, definition_FrameSet)

//-------------------------------------- Definition -----------------------------------------

%include "lsst/meas/multifit/definition/Definition.h"

