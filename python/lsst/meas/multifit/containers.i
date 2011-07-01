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

//-------------------------------------- Set -----------------------------------------

%{
#include "lsst/meas/multifit/containers/Set.h"
%}

namespace lsst { namespace meas { namespace multifit { namespace containers {

template <typename T> class MutableSet {};
template <typename T> class ImmutableSet {};

}}}} // namespace lsst::meas::multifit::containers

%define %DeclareSetBase(T, SET)
%extend SET {
    boost::shared_ptr< T > __getitem__(ID id) {
        return self->get(id);
    }
    int __len__() {
        return self->size();
    }
    PyObject * keys() {
        PyObject * result = PyList_New(0);
        for (SET::iterator i = self->begin(); i != self->end(); ++i) {
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
            return "NAME([\n"  + ",".join(repr(v) for v in self) + "\n])"
    %}
}
%enddef

%define %DeclareImmutableSet(T, NAME)
%template(NAME) lsst::meas::multifit::containers::ImmutableSet< T >;
%DeclareSetBase(T, lsst::meas::multifit::containers::ImmutableSet< T >)
%enddef

%define %DeclareMutableSet(T, NAME)
%template(NAME) lsst::meas::multifit::containers::MutableSet< T >;
%DeclareSetBase(T, lsst::meas::multifit::containers::MutableSet< T >)
%extend lsst::meas::multifit::containers::MutableSet< T > {
    bool insert(boost::shared_ptr<T> const & item) {
        return self->insert(item).second;
    }
    void clear() {
        self->clear();
    }
    void __delitem__(ID id) {
        self->erase(id);
    }
}
%enddef
