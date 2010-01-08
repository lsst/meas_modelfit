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
%typemap(typecheck) ndarray::Array<T,N,C>, ndarray::Array<T,N,C> const * {
    ndarray::PyPyr tmp($input,true);
    $1 = ndarray::PyConverter< ndarray::Array<T,N,C> >::fromPythonStage1(tmp);
    if (!$1) PyErr_Clear();
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

%define %declareEigenMatrix(TYPE)
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
