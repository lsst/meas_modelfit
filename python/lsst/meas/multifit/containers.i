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
#include "lsst/meas/multifit/containers/Array.h"
%}

%define %ExtendArray()
    %extend {
        int __len__() {
            return self->size();
        }
        boost::shared_ptr< T > _get_internal(int n) {
            return boost::shared_ptr< T >(self->begin() + n);
        }
        %pythoncode %{
            def __iter__(self):
                for i in xrange(len(self)):
                    yield self._get_internal(i)
            def __str__(self):
                return "[\n"  + ",".join(str(v) for v in self) + "\n]"
            def __repr__(self):
                return "[\n"  + ",".join(repr(v) for v in self) + "\n]"
            def keys(self):
                return [v.id for v in self]
            def values(self):
                return [v for v in self]
            def __getitem__(self, id):
                for v in self:
                    if v.id == id:
                        return v
                else:
                    raise KeyErrror("Item with ID %s not found" % id)
        %}
    }
%enddef

namespace lsst { namespace meas { namespace multifit { namespace containers {

template <typename T> class MutableSet {};
template <typename T> class ImmutableSet {};

template <typename T, ArrayIndexEnum index>
class Array {
public:
    %ExtendArray()
};

template <typename T, ArrayIndexEnum index>
class ArrayView {
public:
    %ExtendArray()
};

}}}} // namespace lsst::meas::multifit::containers

%define %DeclareSetBase(T, SET, NAME)
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
        def __str__(self):
            return "NAME([\n"  + ",".join(str(v) for v in self) + "\n])"
        def __repr__(self):
            return "NAME([\n"  + ",".join(repr(v) for v in self) + "\n])"
    %}
}
%enddef

%define %DeclareImmutableSet(T, NAME)
%template(NAME) lsst::meas::multifit::containers::ImmutableSet< T >;
%DeclareSetBase(T, lsst::meas::multifit::containers::ImmutableSet< T >, NAME)
%enddef

%define %DeclareMutableSet(T, NAME)
%template(NAME) lsst::meas::multifit::containers::MutableSet< T >;
%DeclareSetBase(T, lsst::meas::multifit::containers::MutableSet< T >, NAME)
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

%define %DeclareArray(T, NAME, INDEX)
%template(NAME) lsst::meas::multifit::containers::Array< T, lsst::meas::multifit::containers::INDEX >;
%enddef

%define %DeclareArrayView(T, NAME, INDEX)
%template(NAME) lsst::meas::multifit::containers::ArrayView< T, lsst::meas::multifit::containers::INDEX >;
%enddef
