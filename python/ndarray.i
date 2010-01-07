%{
#include "ndarray/python.hpp"
%}

%include "exception.i"

%define %ndarrayTypemaps(T, N, C...)
%template(array) ndarray::Array<T,N,C>;
%typemap(out) ndarray::Array<T,N,C> {
    $result = ndarray::PyConverter<ndarray::Array<T,N,C> >::toPython($1);
}
%typemap(typecheck) ndarray::<T,N,C> {
    // 'tmp' is a temporary defined by the 'in' typemap
    tmp = ndarray::PyPtr($input,true);
    $1 = ndarray::PyConverter< ndarray::Array<T,N,C> >::fromPythonStage1(tmp);
    if (!$1) PyErr_Clear();
}
%typemap(in) ndarray::Array<T,N,C> (ndarray::PyPtr tmp) {
    if (!tmp) {  // it should only be set in advance if the typecheck typemap did it
        tmp = ndarray::PyPtr($input,true);
        if (!ndarray::PyConverter< ndarray::Array<T,N,C> >::fromPythonStage1(tmp)) return NULL;
    }
    if (!ndarray::PyConverter< CLASS >::fromPythonStage2(tmp, $1)) return NULL;
}
%enddef

%define %declareArray(T, N, C...)
%template(array)eups 
%ndarrayTypemaps(T,N,C);
%enddef
