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
#pragma SWIG nowarn=401                 // nothin known about base class X
%{
// these sholdn't be necessary, but SWIG fails if they aren't there.
#include "boost/shared_ptr.hpp"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/footprintUtils.h"
#include "lsst/meas/multifit/WindowedFootprint.h"
#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/ModelProjection.h"
#include "lsst/meas/multifit/ModelEvaluator.h"
#include "lsst/meas/multifit/Cache.h"
#include "lsst/meas/multifit/components/Astrometry.h"
#include "lsst/meas/multifit/components/MorphologyProjection.h"
#include "lsst/meas/multifit/components/Morphology.h"
#include "lsst/meas/multifit/components/FixedNonlinearMorphology.h"
#include "lsst/meas/multifit/components/FixedAstrometry.h"
#include "lsst/meas/multifit/components/FourierMorphologyProjection.h"
#include "lsst/meas/multifit/components/PointSourceMorphology.h"
#include "lsst/meas/multifit/components/PointSourceMorphologyProjection.h"
#include "lsst/meas/multifit/components/SersicMorphology.h"
#include "lsst/meas/multifit/components/SersicMorphologyProjection.h"
#include "lsst/meas/multifit/components/ExponentialMorphology.h"
#include "lsst/meas/multifit/components/ExponentialMorphologyProjection.h"
#include "lsst/meas/multifit/ComponentModel.h"
#include "lsst/meas/multifit/ComponentModelProjection.h"
#include "lsst/meas/multifit/FourierModelProjection.h"
#include "lsst/meas/multifit/SingleLinearParameterFitter.h"
#include "lsst/meas/multifit/MinuitFitter.h"
#include "lsst/meas/multifit/ModelFactory.h"
#define NDARRAY_PYTHON_MAIN
#include "ndarray/python.hpp"
#include "ndarray/python/eigen.hpp"
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
%include "lsst/base.h"
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
            version_eups = eups.getSetupVersion("meas_multifit")
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
%import "lsst/afw/coord/coordLib.i"

%include "lsst/afw/image/lsstImageTypes.i"     // Image/Mask types and typedefs

%include "lsst/meas/multifit/core.h"



%define %downcast(BaseType, DerivedType...)
   %extend DerivedType {
       static boost::shared_ptr<DerivedType > swigConvert(
           boost::shared_ptr<BaseType> const & ptr
       ) {
           return boost::dynamic_pointer_cast<DerivedType >(ptr);
       }
   }
%enddef 


%template(ExposureListF) std::list<lsst::afw::image::Exposure<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>::Ptr>;
%template(ExposureListD) std::list<lsst::afw::image::Exposure<double, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>::Ptr>;
%extend lsst::afw::image::Exposure<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel> {
    %pythoncode {
        def makeList(self): 
            return ExposureListF()
    }
}
%extend lsst::afw::image::Exposure<double, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel> {
    %pythoncode {
        def makeList(self):
            return ExposureListD()
    }
}

%declareArray(double, 2, 0);
%declareArray(double, 1, 1);
%declareArray(double const, 1, 1);

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

SWIG_SHARED_PTR(WindowedFootprintPtr, lsst::meas::multifit::WindowedFootprint)
%include "lsst/meas/multifit/WindowedFootprint.h"
%extend lsst::meas::multifit::WindowedFootprint {
    %template(compress) compress<double, double, 0>; 

    %template(expand) expand<double, double, 0>; 
};

SWIG_SHARED_PTR(CahcePtr, lsst::meas::multifit::Cache)
%include "lsst/meas/multifit/Cache.h"

SWIG_SHARED_PTR(ModelPtr, lsst::meas::multifit::Model);
%include "lsst/meas/multifit/Model.h"

SWIG_SHARED_PTR(ModelProjectionPtr, lsst::meas::multifit::ModelProjection)
%include "lsst/meas/multifit/ModelProjection.h"
%extend lsst::meas::multifit::ModelProjection {
    %returnArray(computeModelImage, lsst::meas::multifit::Pixel const, 1, 1);
    %returnArray(computeLinearParameterDerivative, lsst::meas::multifit::Pixel const, 2, 1);
    %returnArray(computeNonlinearParameterDerivative, lsst::meas::multifit::Pixel const, 2, 1);
    //%returnArray(computeWcsParameterDerivative, lsst::meas::multifit::Pixel const, 2, 1);
    //%returnArray(computePsfParameterDerivative, lsst::meas::multifit::Pixel const, 2, 1);
};

SWIG_SHARED_PTR(AstrometryPtr, lsst::meas::multifit::components::Astrometry);
SWIG_SHARED_PTR(MorphologyPtr, lsst::meas::multifit::components::Morphology);
SWIG_SHARED_PTR_DERIVED(FixedAstrometryPtr, lsst::meas::multifit::components::Astrometry,
    lsst::meas::multifit::components::FixedAstrometry);
SWIG_SHARED_PTR_DERIVED(FixedNonlinearMorphologyPtr, lsst::meas::multifit::components::Morphology,
    lsst::meas::multifit::components::FixedNonlinearMorphology);    
SWIG_SHARED_PTR_DERIVED(PointSourceMorphologyPtr, lsst::meas::multifit::components::Morphology,
    lsst::meas::multifit::components::PointSourceMorphology);   
SWIG_SHARED_PTR_DERIVED(SersicMorphologyPtr, lsst::meas::multifit::components::Morphology,
    lsst::meas::multifit::components::SersicMorphology);    
SWIG_SHARED_PTR_DERIVED(ExponentialMorphologyPtr, lsst::meas::multifit::components::Morphology,
    lsst::meas::multifit::components::ExponentialMorphology);    
SWIG_SHARED_PTR(MorphologyProjectionPtr, lsst::meas::multifit::components::MorphologyProjection);
SWIG_SHARED_PTR_DERIVED(FourierMorphologyProjectionPtr,
    lsst::meas::multifit::components::MorphologyProjection,
    lsst::meas::multifit::components::FourierMorphologyProjection);
SWIG_SHARED_PTR_DERIVED(PointSourceMorphologyProjectionPtr,
    lsst::meas::multifit::components::FourierMorphologyProjection,
    lsst::meas::multifit::components::PointSourceMorphologyProjection);
SWIG_SHARED_PTR_DERIVED(SersicMorphologyProjectionPtr,
    lsst::meas::multifit::components::FourierMorphologyProjection,
    lsst::meas::multifit::components::SersicMorphologyProjection);
SWIG_SHARED_PTR_DERIVED(ExponentialMorphologyProjectionPtr,
    lsst::meas::multifit::components::FourierMorphologyProjection,
    lsst::meas::multifit::components::ExponentialMorphologyProjection);
%ignore lsst::meas::multifit::components::PointSourceMorphology::create;
%ignore lsst::meas::multifit::components::SersicMorphology::create;
%ignore lsst::meas::multifit::components::ExponentialMorphology::create;
%ignore lsst::meas::multifit::components::FixedNonlinearMorphology::create;

%include "lsst/meas/multifit/components/Astrometry.h"
%include "lsst/meas/multifit/components/MorphologyProjection.h"
%include "lsst/meas/multifit/components/Morphology.h"
%include "lsst/meas/multifit/components/PointSourceMorphology.h"
%include "lsst/meas/multifit/components/SersicMorphology.h"
%include "lsst/meas/multifit/components/ExponentialMorphology.h"
%include "lsst/meas/multifit/components/FourierMorphologyProjection.h"
%include "lsst/meas/multifit/components/PointSourceMorphologyProjection.h"
%include "lsst/meas/multifit/components/SersicMorphologyProjection.h"
%include "lsst/meas/multifit/components/ExponentialMorphologyProjection.h"

%include "lsst/meas/multifit/components/FixedAstrometry.h"
%include "lsst/meas/multifit/components/FixedNonlinearMorphology.h"

%inline %{
  
    boost::shared_ptr<lsst::meas::multifit::components::PointSourceMorphology> createPointSourceMorphology(
        lsst::meas::multifit::Parameter flux
    ) {
        return lsst::meas::multifit::components::PointSourceMorphology::create(flux);
    }
    
    boost::shared_ptr<lsst::meas::multifit::components::SersicMorphology> createSersicMorphology(
        lsst::meas::multifit::Parameter flux,
        lsst::afw::geom::ellipses::Core const & ellipse,
        lsst::meas::multifit::Parameter sersicIndex
    ) {
        return lsst::meas::multifit::components::SersicMorphology::create(
            flux, ellipse, sersicIndex
        );
    }
    
    boost::shared_ptr<lsst::meas::multifit::components::ExponentialMorphology> createExponentialMorphology(
        lsst::meas::multifit::Parameter flux,
        lsst::afw::geom::ellipses::Core const & ellipse
    ) {
        return lsst::meas::multifit::components::ExponentialMorphology::create(flux, ellipse);
    }

    boost::shared_ptr<lsst::meas::multifit::components::FixedNonlinearMorphology> createFixedNonlinearMorphology(
        lsst::meas::multifit::components::Morphology const & base
    ) {
        return lsst::meas::multifit::components::FixedNonlinearMorphology::create(base);
    }
    
    boost::shared_ptr<lsst::meas::multifit::Model> createSersicModel(
        lsst::meas::multifit::Parameter flux, 
        lsst::afw::geom::Point2D const & centroid,
        lsst::afw::geom::ellipses::Core const & ellipse,
        lsst::meas::multifit::Parameter sersicIndex
    ) {
        return lsst::meas::multifit::ModelFactory::createSersicModel(
            flux, centroid, ellipse, sersicIndex
        );
    }

    boost::shared_ptr<lsst::meas::multifit::Model> createExponentialModel(
        lsst::meas::multifit::Parameter flux, 
        lsst::afw::geom::Point2D const & centroid,
        lsst::afw::geom::ellipses::Core const & ellipse
    ) {
        return lsst::meas::multifit::ModelFactory::createExponentialModel(
            flux, centroid, ellipse
        );
    }

    boost::shared_ptr<lsst::meas::multifit::Model> createPointSourceModel(
        lsst::meas::multifit::Parameter flux, 
        lsst::afw::geom::Point2D const & centroid
    ) {
        return lsst::meas::multifit::ModelFactory::createPointSourceModel(
            flux, centroid
        );
    }
%}

SWIG_SHARED_PTR_DERIVED(ComponentModelPtr, lsst::meas::multifit::Model,     
    lsst::meas::multifit::ComponentModel);
%include "lsst/meas/multifit/ComponentModel.h"

%downcast(lsst::meas::multifit::Model, lsst::meas::multifit::ComponentModel);

SWIG_SHARED_PTR_DERIVED(ComponentModelProjectionPtr, lsst::meas::multifit::ModelProjection,
    lsst::meas::multifit::ComponentModelProjection);
SWIG_SHARED_PTR_DERIVED(FourierModelProjectionPtr, lsst::meas::multifit::ComponentModelProjection,
    lsst::meas::multifit::FourierModelProjection);
%include "lsst/meas/multifit/ComponentModelProjection.h"    
%include "lsst/meas/multifit/FourierModelProjection.h"

%downcast(lsst::meas::multifit::ModelProjection, lsst::meas::multifit::FourierModelProjection);
%downcast(lsst::meas::multifit::ModelProjection, lsst::meas::multifit::ComponentModelProjection);
%downcast(lsst::meas::multifit::ComponentModelProjection, lsst::meas::multifit::FourierModelProjection);

SWIG_SHARED_PTR(ModelEvaluatorPtr, lsst::meas::multifit::ModelEvaluator);
%nodefaultctor lsst::meas::multifit::ModelEvaluator;
%include "lsst/meas/multifit/ModelEvaluator.h"
%extend lsst::meas::multifit::ModelEvaluator {
    %returnArray(getDataVector, lsst::meas::multifit::Pixel const, 1, 1);
    
    %template(setExposureList) setExposureList<double, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>; 
    %template(setExposureList) setExposureList<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>;
};

SWIG_SHARED_PTR(SimpleResultPtr, lsst::meas::multifit::SimpleFitResult);
%include "lsst/meas/multifit/SingleLinearParameterFitter.h"


%include "lsst/meas/multifit/MinuitFitter.h"
