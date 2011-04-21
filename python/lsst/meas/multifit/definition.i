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

//----------------------------- ParameterComponents -----------------------------------

%{
#include "lsst/meas/multifit/definition/parameters.h"
%}

SWIG_SHARED_PTR(definition_PositionComponentPtr, lsst::meas::multifit::definition::ParameterComponent<POSITION>);
SWIG_SHARED_PTR(definition_RadiusComponentPtr, lsst::meas::multifit::definition::ParameterComponent<RADIUS>);
SWIG_SHARED_PTR(definition_EllipticityComponentPtr, lsst::meas::multifit::definition::ParameterComponent<ELLIPTICITY>);
%AddStreamRepr(lsst::meas::multifit::definition::ParameterComponent<POSITION>)
%AddStreamRepr(lsst::meas::multifit::definition::ParameterComponent<RADIUS>)
%AddStreamRepr(lsst::meas::multifit::definition::ParameterComponent<ELLIPTICITY>)

%include "lsst/meas/multifit/definition/parameters.h"

%define %ParameterComponent_POSTINCLUDE(NAME, ENUM)
%template(definition_##NAME##Component)
lsst::meas::multifit::definition::ParameterComponent<lsst::meas::multifit::ENUM>;
%PointerEQ(lsst::meas::multifit::definition::ParameterComponent<lsst::meas::multifit::ENUM>)
%enddef

%ParameterComponent_POSTINCLUDE(Position, POSITION);
%ParameterComponent_POSTINCLUDE(Radius, RADIUS);
%ParameterComponent_POSTINCLUDE(Ellipticity, ELLIPTICITY);

//------------------------------------- Object ---------------------------------------

%{
#include "lsst/meas/multifit/definition/Object.h"
%}

SWIG_SHARED_PTR(definition_ObjectPtr, lsst::meas::multifit::definition::Object);

%rename(definition_Object) lsst::meas::multifit::definition::Object;

%include "lsst/meas/multifit/definition/Object.h"

%PointerEQ(lsst::meas::multifit::definition::Object)
%AddStreamRepr(lsst::meas::multifit::definition::Object)

//------------------------------------- Frame ---------------------------------------

%{
#include "lsst/meas/multifit/definition/Frame.h"
%}

SWIG_SHARED_PTR(definition_FramePtr, lsst::meas::multifit::definition::Frame);

%rename(definition_Frame) lsst::meas::multifit::definition::Frame;

%include "lsst/meas/multifit/definition/Frame.h"

%PointerEQ(lsst::meas::multifit::definition::Frame)
%AddStreamRepr(lsst::meas::multifit::definition::Frame)

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

%DeclareSet(Object)
%DeclareSet(Frame)


//-------------------------------------- Definition -----------------------------------------

%include "lsst/meas/multifit/Definition.h"
