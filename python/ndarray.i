%{
#include "ndarray/python.hpp"
%}

%include "exception.i"

%define ndarrayTypemaps(CLASS)
%typemap(out) CLASS {
    $result = ndarray::PyConverter< CLASS >::toPython($1);
}
%typemap(typecheck) CLASS {
    // 'tmp' is a temporary defined by the 'in' typemap
    tmp = ndarray::PyPtr($input,true);
    $1 = ndarray::PyConverter< CLASS >::fromPythonStage1(tmp);
    if (!$1) PyErr_Clear();
}
%typemap(in) CLASS (ndarray::PyPtr tmp) {
    if (!tmp) {  // it should only be set in advance if the typecheck typemap did it
        tmp = ndarray::PyPtr($input,true);
        if (!ndarray::PyConverter< CLASS >::fromPythonStage1(tmp)) return NULL;
    }
    if (!ndarray::PyConverter< CLASS >::fromPythonStage2(tmp, $1)) return NULL;
}
%enddef

%define declareArray(T, N, C)
%ndarrayTypemaps(ndarray::Array<T,N,C>)
%enddef
