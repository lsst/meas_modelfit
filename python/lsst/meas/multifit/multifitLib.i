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
#pragma SWIG nowarn=401					// nothin known about base class X
%{
// these sholdn't be necessary, but SWIG fails if they aren't there.
#include "lsst/afw/detection.h" 
#include "lsst/meas/algorithms/Centroid.h"
#include "lsst/meas/algorithms/Photometry.h"
#include "lsst/meas/algorithms/Shape.h"

#include "lsst/afw/detection/Footprint.h"
#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/footprintUtils.h"
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
#define NDARRAY_PYTHON_MAIN
%}

%inline %{
namespace boost { }
namespace lsst { namespace meas { namespace multifit { namespace components {} } } }    
%}

%ignore boost::noncopyable;
namespace boost {
    class noncopyable {};
}

%include "std_list.i"
%include "ndarray.i"

%{
#include "lsst/afw/numpyTypemaps.h"
%}
%init %{
    import_array();
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
%import "lsst/afw/geom/geomLib.i"
%import "lsst/afw/image/imageLib.i"
%import "lsst/afw/detection/detectionLib.i"
%import "lsst/afw/math/mathLib.i"


%import "lsst/meas/algorithms/algorithmsLib.i"

%include "lsst/meas/multifit/core.h"

%define %downcast(BaseType, DerivedType...)
   %extend DerivedType {
       static boost::shared_ptr<DerivedType > swigConvert(boost::shared_ptr<BaseType> const & ptr) {
           return boost::dynamic_pointer_cast<DerivedType >(ptr);
       }
   }
%enddef 

%declareArray(lsst::meas::multifit::Pixel const, 1, 1);
%declareArray(lsst::meas::multifit::Pixel const, 2, 1);
%declareArray(lsst::meas::multifit::Pixel, 1, 1);
%declareArray(lsst::meas::multifit::Pixel, 2, 1);
%declareEigenMatrix(lsst::meas::multifit::ParameterVector);

%include "lsst/meas/multifit/footprintUtils.h"
%template(compressImageF) lsst::meas::multifit::compressImage<float, 
    lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>;
%template(compressImageD) lsst::meas::multifit::compressImage<double, 
    lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>;
%template(expandImageF) lsst::meas::multifit::expandImage<float, 
    lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>;
%template(expandImageD) lsst::meas::multifit::expandImage<double, 
    lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>;

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

%downcast(lsst::meas::multifit::ModelFactory, lsst::meas::multifit::ComponentModelFactory);
%downcast(lsst::meas::multifit::ModelFactory, lsst::meas::multifit::PointSourceModelFactory);
%downcast(lsst::meas::multifit::ComponentModelFactory, lsst::meas::multifit::PointSourceModelFactory);


SWIG_SHARED_PTR(ModelProjectionPtr, lsst::meas::multifit::ModelProjection)
%include "lsst/meas/multifit/ModelProjection.h"
%extend lsst::meas::multifit::ModelProjection {
    %returnArray(computeModelImage, lsst::meas::multifit::Pixel const, 1, 1);
    %returnArray(computeLinearParameterDerivative, lsst::meas::multifit::Pixel const, 2, 1);
    %returnArray(computeNonlinearParameterDerivative, lsst::meas::multifit::Pixel const, 2, 1);
    %returnArray(computeWcsParameterDerivative, lsst::meas::multifit::Pixel const, 2, 1);
    %returnArray(computePsfParameterDerivative, lsst::meas::multifit::Pixel const, 2, 1);
};


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

%downcast(lsst::meas::multifit::Model, lsst::meas::multifit::ComponentModel);

SWIG_SHARED_PTR_DERIVED(ComponentModelProjectionPtr, lsst::meas::multifit::ModelProjection,
    lsst::meas::multifit::ComponentModelProjection)
SWIG_SHARED_PTR_DERIVED(FourierModelProjectionPtr, lsst::meas::multifit::ComponentModelProjection,
    lsst::meas::multifit::FourierModelProjection)
%include "lsst/meas/multifit/ComponentModelProjection.h"    
%include "lsst/meas/multifit/FourierModelProjection.h"

%downcast(lsst::meas::multifit::ModelProjection, lsst::meas::multifit::FourierModelProjection);
%downcast(lsst::meas::multifit::ModelProjection, lsst::meas::multifit::ComponentModelProjection);
%downcast(lsst::meas::multifit::ComponentModelProjection, lsst::meas::multifit::FourierModelProjection);

SWIG_SHARED_PTR_DERIVED(
    CharacterizedExposureFPtr, 
    lsst::afw::image::Exposure<float>,
    lsst::meas::multifit::CharacterizedExposure<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>
);                      
SWIG_SHARED_PTR_DERIVED(
    CharacterizedExposureDPtr, 
    lsst::afw::image::Exposure<double>,
    lsst::meas::multifit::CharacterizedExposure<double, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>
); 

%include "lsst/meas/multifit/CharacterizedExposure.h"


SWIG_SHARED_PTR(ModelEvaluatorPtr, lsst::meas::multifit::ModelEvaluator)
%nodefaultctor lsst::meas::multifit::ModelEvaluator;
%include "lsst/meas/multifit/ModelEvaluator.h"
%extend lsst::meas::multifit::ModelEvaluator {
    %returnArray(getImageVector, lsst::meas::multifit::Pixel const, 1, 1);
    %returnArray(getVarianceVector, lsst::meas::multifit::Pixel const, 1, 1);
    %returnArray(computeModelImage, lsst::meas::multifit::Pixel const, 1, 1);
    %returnArray(computeLinearParameterDerivative, lsst::meas::multifit::Pixel const, 2, 2);
    %returnArray(computeNonlinearParameterDerivative, lsst::meas::multifit::Pixel const, 2, 2);
    
    %template(ModelEvaluator) ModelEvaluator<double, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>; 
    %template(ModelEvaluator) ModelEvaluator<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>; 
    %template(setExposureList) setExposureList<double, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>; 
    %template(setExposureList) setExposureList<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>;
};

SWIG_SHARED_PTR(SimpleResultPtr, lsst::meas::multifit::SimpleFitResult)
%include "lsst/meas/multifit/SingleLinearParameterFitter.h"

%template(CharacterizedExposureF) lsst::meas::multifit::CharacterizedExposure<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>;
%template(CharacterizedExposureD) lsst::meas::multifit::CharacterizedExposure<double, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>;
%template(CharacterizedExposureListF) std::list<lsst::meas::multifit::CharacterizedExposure<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>::Ptr>;
%template(CharacterizedExposureListD) std::list<lsst::meas::multifit::CharacterizedExposure<double, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>::Ptr>;
