%{
#include "ndarray/python.hpp"
%}

%define %declareArray(T,N,C)
// BEGIN BROKEN =================================================
//%typemap(out, optimal="1") ndarray::Array<T,N,C> {
//    $result = ndarray::PyConverter< ndarray::Array<T,N,C> >::toPython($1);
//}
// END BROKEN ===================================================
%typemap(typecheck) ndarray::Array<T,N,C>, ndarray::Array<T,N,C> const & {
    // 'tmp' is a temporary defined by the 'in' typemap
    tmp = ndarray::PyPtr($input,true);
    $1 = ndarray::PyConverter< ndarray::Array<T,N,C> >::fromPythonStage1(tmp);
    if (!$1) PyErr_Clear();
}
%typemap(in) ndarray::Array<T,N,C>, ndarray::Array<T,N,C> const & (ndarray::PyPtr tmp) {
    if (!tmp) {  // it should only be set in advance if the typecheck typemap did it
        tmp = ndarray::PyPtr($input,true);
        if (!ndarray::PyConverter< ndarray::Array<T,N,C> >::fromPythonStage1(tmp)) return NULL;
    }
    if (!ndarray::PyConverter< ndarray::Array<T,N,C> >::fromPythonStage2(tmp, $1)) return NULL;
}
%enddef

%define %returnArray(METHOD, T, N, C)
PyObject * METHOD() {
    ndarray::Array<T,N,C> r(self->METHOD());
    return ndarray::PyConverter< ndarray::Array<T,N,C> >::toPython(r);
}
%enddef
