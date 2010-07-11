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
 
%{
#include "ndarray/python.hpp"
#include "ndarray/python/eigen.hpp"
#include <boost/scoped_ptr.hpp>
%}

%define %declareArray(T,N,C)
/* BROKEN ======================================================================
// SWIG assumes default construction + assignment are available and are
// equivalent to copy construction; but that's not true for ndarray.
%typemap(out, optimal="1") ndarray::Array<T,N,C> {
    $result = ndarray::PyConverter< ndarray::Array<T,N,C> >::toPython($1);
}
============================================================================= */
%typemap(typecheck) ndarray::Array<T,N,C>, ndarray::Array<T,N,C> const *, ndarray::Array<T,N,C> const & {
    ndarray::PyPtr tmp($input,true);
    $1 = ndarray::PyConverter< ndarray::Array<T,N,C> >::fromPythonStage1(tmp);
    //if (!$1) PyErr_Clear();
    return NULL;
}
%typemap(in) ndarray::Array<T,N,C> const & (ndarray::Array<T,N,C> val) {
    ndarray::PyPtr tmp($input,true);
    if (!ndarray::PyConverter< ndarray::Array<T,N,C> >::fromPythonStage1(tmp)) return NULL;
    if (!ndarray::PyConverter< ndarray::Array<T,N,C> >::fromPythonStage2(tmp, val)) return NULL;
    $1 = &val;
}
%typemap(in) ndarray::Array<T,N,C> {
    ndarray::PyPtr tmp($input,true);
    if (!ndarray::PyConverter< ndarray::Array<T,N,C> >::fromPythonStage1(tmp)) return NULL;
    if (!ndarray::PyConverter< ndarray::Array<T,N,C> >::fromPythonStage2(tmp, $1)) return NULL;
}
%enddef

// Partial workaround for the broken output typemap above; only works for methods with no arguments.
%define %returnArray(METHOD, T, N, C)
PyObject * METHOD() {
    ndarray::Array<T,N,C> r(self->METHOD());
    return ndarray::PyConverter< ndarray::Array<T,N,C> >::toPython(r);
}
%enddef

%define %eigenMatrix(TYPE)
%typemap(out) TYPE {
    $result = ndarray::PyConverter< TYPE >::toPython($1);
}
%typemap(out) TYPE const & {
    $result = ndarray::PyConverter< TYPE >::toPython(*$1);
}
%typemap(typecheck) TYPE, TYPE const * {
    ndarray::PyPtr tmp($input,true);
    $1 = ndarray::PyConverter< TYPE >::fromPythonStage1(tmp);
    if (!$1) PyErr_Clear();
}
%typemap(in) TYPE const & (TYPE val) {
    ndarray::PyPtr tmp($input,true);
    if (!ndarray::PyConverter< TYPE >::fromPythonStage1(tmp)) return NULL;
    if (!ndarray::PyConverter< TYPE >::fromPythonStage2(tmp, val)) return NULL;
    $1 = &val;
}
%enddef
