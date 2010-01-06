// -*- lsst-c++ -*-
%define multifitLib_DOCSTRING
"
Basic routines to talk to lsst::meas::multifit classes
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.meas.multifit", docstring=multifitLib_DOCSTRING) multifitLib

// Suppress swig complaints
#pragma SWIG nowarn=314                 // print is a python keyword (--> _print)
#pragma SWIG nowarn=362                 // operator=  ignored

%{
#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/ModelProjection.h"
#include "lsst/meas/multifit/ModelFactory.h"
#include "lsst/meas/multifit/ModelEvaluator.h"
%}

%init %{
%}
/************************************************************************************************************/

%include "lsst/p_lsstSwig.i"
%include "std_complex.i"
	
%lsst_exceptions();

%pythoncode %{
import lsst.utils

def version(HeadURL = r"$HeadURL: svn+ssh://svn.lsstcorp.org/DMS/meas/multifit/trunk/python/lsst/meas/multifit/multifitLib.i $"):
    """Return a version given a HeadURL string. If a different version is setup, return that too"""

    version_svn = lsst.utils.guessSvnVersion(HeadURL)

    try:
        import eups
    except ImportError:
        return version_svn
    else:
        try:
            version_eups = eups.getSetupVersion("meas")
        except AttributeError:
            return version_svn

    if version_eups == version_svn:
        return version_svn
    else:
        return "%s (setup: %s)" % (version_svn, version_eups)
%}

%import "lsst/afw/image/imageLib.i"
%import "lsst/afw/detection/detectionLib.i"
%import "lsst/afw/math/mathLib.i"
%import "lsst/afw/geom/geomLib.i"


%import "lsst/meas/algorithms/algorithmsLib.i"

SWIG_SHARED_PTR(ModelPtr, lsst::meas::multifit::Model);
%include "lsst/meas/multifit/Model.h"

SWIG_SHARED_PTR(ModelProjectionPtr, lsst::meas::multifit::ModelProjection);
%include "lsst/meas/multifit/ModelProjection.h"

SWIG_SHARED_PTR(ModelFactoryPtr, lsst::meas::multifit::ModelFactory);
%include "lsst/meas/multifit/ModelFactory.h"

SWIG_SHARED_PTR(ModelEvalutatorPtr, lsst::meas::multifit::ModelEvaluator);
%include "lsst/meas/multifit/ModelEvaluator.h"
