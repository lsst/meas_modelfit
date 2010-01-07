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
#include "lsst/meas/algorithms/Centroid.h"
#include "lsst/meas/algorithms/Photometry.h"
#include "lsst/meas/algorithms/Shape.h"
#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/ModelProjection.h"
#include "lsst/meas/multifit/ModelEvaluator.h"
#include "lsst/meas/multifit/ModelFactory.h"
#include "lsst/meas/multifit/components/Astrometry.h"
#include "lsst/meas/multifit/components/MorphologyProjection.h"
#include "lsst/meas/multifit/components/Morphology.h"
#include "lsst/meas/multifit/components/FourierMorphologyProjection.h"
#include "lsst/meas/multifit/components/PointSourceMorphology.h"
#include "lsst/meas/multifit/components/PointSourceMorphologyProjection.h"
#include "lsst/meas/multifit/ComponentModel.h"
#include "lsst/meas/multifit/ComponentModelProjection.h"
#include "lsst/meas/multifit/ComponentModelFactory.h"
#include "lsst/meas/multifit/FourierModelProjection.h"
#include "lsst/meas/multifit/PointSourceModelFactory.h"
#include "lsst/meas/multifit/SingleLinearParameterFitter.h"
%}

%inline %{
namespace boost { }
namespace lsst { namespace meas { namespace multifit { namespace components {} } } }    
%}

%ignore boost::noncopyable;
namespace boost {
    class noncopyable {};
}

%init %{
%}
/******************************************************************************/

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

%import "lsst/daf/base/baseLib.i"
%import "lsst/afw/image/imageLib.i"
%import "lsst/afw/detection/detectionLib.i"
%import "lsst/afw/math/mathLib.i"
%import "lsst/afw/geom/geomLib.i"

%import "lsst/meas/algorithms/algorithmsLib.i"

%include "lsst/meas/multifit/core.h"

SWIG_SHARED_PTR(ModelPtr, lsst::meas::multifit::Model)   
%include "lsst/meas/multifit/Model.h"

SWIG_SHARED_PTR(ModelFactoryPtr, lsst::meas::multifit::ModelFactory)
SWIG_SHARED_PTR_DERIVED(ComponentModelFactoryPtr, lsst::meas::multifit::ModelFactory,
    lsst::meas::multifit::ComponentModelFactory)
SWIG_SHARED_PTR_DERIVED(PointSourceModelFactoryPtr, lsst::meas::multifit::ComponentModelFactory,
    lsst::meas::multifit::PointSourceModelFactory)

%include "lsst/meas/multifit/ModelFactory.h"
%include "lsst/meas/multifit/ComponentModelFactory.h"
%include "lsst/meas/multifit/PointSourceModelFactory.h"


SWIG_SHARED_PTR(ModelProjectionPtr, lsst::meas::multifit::ModelProjection)
%ignore lsst::meas::multifit::ModelProjection::computeModelImage;
%ignore lsst::meas::multifit::ModelProjection::computeLinearParameterDerivative;
%ignore lsst::meas::multifit::ModelProjection::computeNonlinearParameterDerivative;
%ignore lsst::meas::multifit::ModelProjection::computeWcsParameterDerivative;
%ignore lsst::meas::multifit::ModelProjection::computePsfParameterDerivative;
%include "lsst/meas/multifit/ModelProjection.h"

SWIG_SHARED_PTR(AstrometryPtr, lsst::meas::multifit::components::Astrometry)
SWIG_SHARED_PTR(MorphologyPtr, lsst::meas::multifit::components::Morphology)
SWIG_SHARED_PTR_DERIVED(PointSourceMorphologyPtr, lsst::meas::multifit::components::Morphology,
    lsst::meas::multifit::components::PointSourceMorphology)    
SWIG_SHARED_PTR(MorphologyProjectionPtr, lsst::meas::multifit::components::MorphologyProjection)
SWIG_SHARED_PTR_DERIVED(FourierMorphologyProjectionPtr,
    lsst::meas::multifit::components::MorphologyProjection,
    lsst::meas::multifit::components::FourierMorphologyProjection)
SWIG_SHARED_PTR_DERIVED(PointSourceMorphologyProjectionPtr,
    lsst::meas::multifit::components::FourierMorphologyProjection,
    lsst::meas::multifit::components::PointSourceMorphologyProjection)
%include "lsst/meas/multifit/components/Astrometry.h"
%include "lsst/meas/multifit/components/MorphologyProjection.h"
%include "lsst/meas/multifit/components/Morphology.h"
%include "lsst/meas/multifit/components/FourierMorphologyProjection.h"
%include "lsst/meas/multifit/components/PointSourceMorphologyProjection.h"
%include "lsst/meas/multifit/components/PointSourceMorphology.h"

SWIG_SHARED_PTR_DERIVED(ComponentModelPtr, lsst::meas::multifit::Model, 	
    lsst::meas::multifit::ComponentModel) 
%include "lsst/meas/multifit/ComponentModel.h"

SWIG_SHARED_PTR_DERIVED(ComponentModelProjectionPtr, lsst::meas::multifit::ModelProjection,
    lsst::meas::multifit::ComponentModelProjection)
SWIG_SHARED_PTR_DERIVED(FourierModelProjectionPtr, lsst::meas::multifit::ComponentModelProjection,
    lsst::meas::multifit::FourierModelProjection)
%include "lsst/meas/multifit/ComponentModelProjection.h"    
%include "lsst/meas/multifit/FourierModelProjection.h"


SWIG_SHARED_PTR(ModelEvaluatorPtr, lsst::meas::multifit::ModelEvaluator)
%ignore lsst::meas::multifit::ModelEvaluator::getImageVector;
%ignore lsst::meas::multifit::ModelEvaluator::getVarianceVector;
%ignore lsst::meas::multifit::ModelEvaluator::computeModelImage;
%ignore lsst::meas::multifit::ModelEvaluator::computeLinearParameterDerivative;
%ignore lsst::meas::multifit::ModelEvaluator::computeNonlinearParameterDerivative;

%include "lsst/meas/multifit/ModelEvaluator.h"

SWIG_SHARED_PTR(SimpleResultPtr, lsst::meas::multifit::SimpleFitResult)
%include "lsst/meas/multifit/SingleLinearParameterFitter.h"

