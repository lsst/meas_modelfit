// -*- lsst-c++ -*-
%define meas_multifit_DOCSTRING
"
Access to the classes from the meas_multifit library
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.meas.multifit", docstring=meas_multifit_DOCSTRING) multifitLib

%include "lsst/p_lsstSwig.i"

%lsst_exceptions()
%import "lsst/pex/exceptions/exceptionsLib.i"

