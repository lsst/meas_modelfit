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

//----------------------------- ParameterElements -----------------------------------

%{
#include "lsst/meas/multifit/definition/elements.h"
%}

SWIG_SHARED_PTR(definition_PositionElementPtr, lsst::meas::multifit::definition::ParameterElement<lsst::meas::multifit::POSITION>);
SWIG_SHARED_PTR(definition_RadiusElementPtr, lsst::meas::multifit::definition::ParameterElement<lsst::meas::multifit::RADIUS>);
SWIG_SHARED_PTR(definition_EllipticityElementPtr, lsst::meas::multifit::definition::ParameterElement<lsst::meas::multifit::ELLIPTICITY>);
SWIG_SHARED_PTR(definition_FluxElementPtr, lsst::meas::multifit::definition::FluxElement);
%rename(definition_FluxElement) lsst::meas::multifit::definition::FluxElement;

%define %DeclareDefinitionParameterElement(TITLE, LOWER, UPPER, CONSTRAINT)
%template(definition_##TITLE##Element)
lsst::meas::multifit::definition::ParameterElement<lsst::meas::multifit::UPPER>;
%PointerEQ(lsst::meas::multifit::definition::ParameterElement<lsst::meas::multifit::UPPER>)
%AddStreamRepr(lsst::meas::multifit::definition::ParameterElement<lsst::meas::multifit::UPPER>)
%extend lsst::meas::multifit::definition::ParameterElement<lsst::meas::multifit::UPPER> {
    static boost::shared_ptr< 
        lsst::meas::multifit::definition::ParameterElement< lsst::meas::multifit::UPPER >
    > make(lsst::meas::multifit::TITLE const & value, bool active=true) {
        return lsst::meas::multifit::definition::ParameterElement< lsst::meas::multifit::UPPER >::make(
            value, active
        );
    }
    static boost::shared_ptr< 
        lsst::meas::multifit::definition::ParameterElement< lsst::meas::multifit::UPPER >
    > make(lsst::meas::multifit::TITLE const & value, lsst::meas::multifit::CONSTRAINT const & bounds, 
           bool active=true
    ) {
        return lsst::meas::multifit::definition::ParameterElement< lsst::meas::multifit::UPPER >::make(
            value, bounds, active
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
    boost::shared_ptr< lsst::meas::multifit::definition::ParameterElement< lsst::meas::multifit::UPPER > > get##TITLE##Element() {
        return self->get##TITLE##Element();
    }
    void set##TITLE##Element(boost::shared_ptr< lsst::meas::multifit::definition::ParameterElement< lsst::meas::multifit::UPPER > > value) {
        self->get##TITLE##Element() = value;
    }
}
%enddef
%rename(detail_CircleConstraint) lsst::meas::multifit::detail::CircleConstraint;
%rename(detail_MinMaxConstraint) lsst::meas::multifit::detail::MinMaxConstraint;

%include "lsst/meas/multifit/definition/elements.h"

//------------------------------------- ObjectComponent ---------------------------------------

%{
#include "lsst/meas/multifit/definition/ObjectComponent.h"
%}

SWIG_SHARED_PTR(detail_ObjectComponentBasePtr, lsst::meas::multifit::detail::ObjectComponentBase);
SWIG_SHARED_PTR_DERIVED(definition_ObjectComponentPtr, lsst::meas::multifit::detail::ObjectComponentBase, lsst::meas::multifit::definition::ObjectComponent);

%rename(definition_ObjectComponent) lsst::meas::multifit::definition::ObjectComponent;

%include "lsst/meas/multifit/definition/ObjectComponent.h"

%DeclareDefinitionParameterElement(Position, position, POSITION, detail::CircleConstraint);
%DeclareDefinitionParameterElement(Radius, radius, RADIUS, detail::MinMaxConstraint);
%DeclareDefinitionParameterElement(Ellipticity, ellipticity, ELLIPTICITY, detail::CircleConstraint);

%PointerEQ(lsst::meas::multifit::definition::FluxElement)
%AddStreamRepr(lsst::meas::multifit::definition::FluxElement)
%extend lsst::meas::multifit::definition::FluxElement {
    static boost::shared_ptr<lsst::meas::multifit::definition::FluxElement> make() {
        return lsst::meas::multifit::definition::FluxElement::make();
    }
    static boost::shared_ptr<lsst::meas::multifit::definition::FluxElement> make(double value) {
        return lsst::meas::multifit::definition::FluxElement::make(value);
    }
    static boost::shared_ptr<lsst::meas::multifit::definition::FluxElement> make(
        double value, lsst::meas::multifit::detail::MinMaxConstraint const & bounds
    ) {
        return lsst::meas::multifit::definition::FluxElement::make(value, bounds);
    }
    double getValue() { return self->getValue(); }
    void setValue(double value) { self->getValue() = value; }

    lsst::meas::multifit::detail::MinMaxConstraint & getBounds() { return self->getBounds(); }
    void setBounds(lsst::meas::multifit::detail::MinMaxConstraint const & bounds) {
        self->getBounds() = bounds;
    }
}
%extend lsst::meas::multifit::definition::ObjectComponent {
    boost::shared_ptr<lsst::meas::multifit::definition::FluxElement> getFluxElement() {
        return self->getFluxElement();
    }
    void setFluxElement(boost::shared_ptr<lsst::meas::multifit::definition::FluxElement> value) {
        self->getFluxElement() = value;
    }
}

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

%{
#include "lsst/meas/multifit/definition/Set.h"
%}

namespace lsst { namespace meas { namespace multifit { namespace definition {

template <typename Value_> class Set {};

}}}} // namespace lsst::meas::multifit::definition

%define %DeclareSet(NAME)
%template(definition_##NAME##Set)
lsst::meas::multifit::definition::Set<lsst::meas::multifit::definition::NAME>;
%extend lsst::meas::multifit::definition::Set<lsst::meas::multifit::definition::NAME> {
    boost::shared_ptr< NAME > __getitem__(ID id) {
        lsst::meas::multifit::definition::Set<lsst::meas::multifit::definition::NAME>::iterator i 
            = self->find(id);
        boost::shared_ptr< lsst::meas::multifit::definition::NAME > result;
        if (i != self->end()) {
            result = i;
        }
        return result;
    }
    bool insert(boost::shared_ptr< lsst::meas::multifit::definition::NAME > const & item) {
        return self->insert(item).second;
    }
    int __len__() {
        return self->size();
    }
    void clear() {
        self->clear();
    }

    PyObject * keys() {
        PyObject * result = PyList_New(0);
        typedef lsst::meas::multifit::definition::Set<lsst::meas::multifit::definition::NAME> Self;
        for (Self::iterator i = self->begin(); i != self->end(); ++i) {
            PyObject * element = PyLong_FromLongLong(i->id);
            if (element == 0) {
                Py_DECREF(result);
                return 0;
            }
            if (PyList_Append(result, element) < 0) {
                Py_DECREF(result);
                Py_DECREF(element);
                return 0;
            }
            Py_DECREF(element);
        }
        return result;
    }

    %pythoncode %{
        def values(self):
            return [self[id] for id in self.keys()]
        def __iter__(self):
            for id in self.keys():
                yield self[id]
        def __repr__(self):
            return "NAME""Set([\n"  + ",".join(repr(v) for v in self) + "\n])"
    %}
}
%enddef

%DeclareSet(ObjectComponent)
%DeclareSet(Frame)


//-------------------------------------- Definition -----------------------------------------

%include "lsst/meas/multifit/definition/Definition.h"
