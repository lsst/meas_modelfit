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
#include "lsst/pex/logging.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/detection/AperturePhotometry.h"
#include "lsst/afw/cameraGeom.h"
#include "lsst/meas/multifit/constants.h"
#include "lsst/meas/multifit/BaseEvaluator.h"
#include "lsst/meas/multifit/Evaluator.h"
#include "lsst/meas/multifit/Evaluation.h"
#include "lsst/meas/multifit/ModelBasis.h"
#include "lsst/meas/multifit/ShapeletModelBasis.h"
#include "lsst/meas/multifit/CompoundShapeletModelBasis.h"    
#include "lsst/ndarray/eigen.h"
#include <Eigen/Core>
#define PY_ARRAY_UNIQUE_SYMBOL LSST_MEAS_MULTIFIT_NUMPY_ARRAY_API
#include "numpy/arrayobject.h"
#include "lsst/ndarray/python.h"
#include "lsst/ndarray/python/eigen.h"
%}

/******************************************************************************/
%init %{
    import_array();
%}

%include "lsst/p_lsstSwig.i"
%include "lsst/base.h"
%include "std_complex.i"

%lsst_exceptions();

%pythoncode %{
import lsst.utils

def version(HeadURL = r"$HeadURL$"):
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

def makeSourceMeasurement(**kw):
    policy = lsst.pex.policy.Policy()
    algorithm = "SHAPELET_MODEL"
    policy.add("%s.enabled" % algorithm, True)
    for k in kw:
        policy.add("%s.%s" % (algorithm, k), kw[k])
    options = lsst.meas.multifit.SourceMeasurement.readPolicy(policy.get(algorithm))
    measurement = lsst.meas.multifit.SourceMeasurement(options)
    return measurement, policy

%}

%include "lsst/ndarray/ndarray.i"
%import "lsst/afw/geom/geomLib.i"
%import "lsst/afw/geom/ellipses/ellipsesLib.i"
%import "lsst/afw/detection/detectionLib.i"
%import "lsst/afw/math/mathLib.i"
%import "lsst/afw/math/shapelets/shapeletsLib.i"
%import "lsst/afw/image/imageLib.i"

/*****************************************************************************/
%declareNumPyConverters(lsst::ndarray::Array<lsst::meas::multifit::Pixel, 2, 1>);
%declareNumPyConverters(lsst::ndarray::Array<lsst::meas::multifit::Pixel const, 2, 1>);
%declareNumPyConverters(lsst::ndarray::Array<lsst::meas::multifit::Pixel, 2, 2>);
%declareNumPyConverters(lsst::ndarray::Array<lsst::meas::multifit::Pixel const, 2, 2>);
%declareNumPyConverters(lsst::ndarray::Array<lsst::meas::multifit::Pixel, 3, 3>);
%declareNumPyConverters(lsst::ndarray::Array<lsst::meas::multifit::Pixel const, 3, 3>);
%declareNumPyConverters(lsst::ndarray::Array<lsst::meas::multifit::Pixel, 1, 1>);
%declareNumPyConverters(lsst::ndarray::Array<lsst::meas::multifit::Pixel const, 1, 1>);
%declareNumPyConverters(lsst::ndarray::Array<double const, 1, 1>);
%declareNumPyConverters(lsst::ndarray::Array<double, 1, 1>);
%declareNumPyConverters(lsst::ndarray::Array<double, 2, 2>);
%declareNumPyConverters(lsst::ndarray::Array<double const, 2, 2>);
%declareNumPyConverters(lsst::ndarray::Array<double const, 3, 3>);
%declareNumPyConverters(Eigen::VectorXd);

%include "lsst/meas/multifit/constants.h"

%shared_ptr(lsst::meas::multifit::ModelBasis);
%shared_ptr(lsst::meas::multifit::ShapeletModelBasis);
%shared_ptr(lsst::meas::multifit::CompoundShapeletModelBasis);

%nodefaultctor lsst::meas::multifit::ModelBasis;
%nodefaultctor lsst::meas::multifit::ShapeletModelBasis;
%nodefaultctor lsst::meas::multifit::CompoundShapeletModelBasis;


%extend lsst::meas::multifit::CompoundShapeletModelBasis {
    %feature("shadow") _getMapping %{
        def getMapping(self):
            return $action(self)
    %}
    %feature("shadow") _extractComponents %{
        def extractComponents(self):
            return $action(self)
    %}


    lsst::ndarray::Array<lsst::meas::multifit::Pixel const, 2, 1> _getMapping() const {
        return self->getMapping();
    }
    lsst::meas::multifit::CompoundShapeletModelBasis::ComponentVector _extractComponents() const {
        return self->extractComponents();
    }
};

%shared_ptr(lsst::meas::multifit::ProfileFunction);

%include "lsst/meas/multifit/ModelBasis.h"
%include "lsst/meas/multifit/ShapeletModelBasis.h"
%include "lsst/meas/multifit/CompoundShapeletModelBasis.h"

%template(CompoundShapelet_ComponentVector) std::vector<boost::shared_ptr<lsst::meas::multifit::ShapeletModelBasis> >;

%include "definition.i"
%include "grid.i"

%shared_ptr(lsst::meas::multifit::BaseEvaluator);

%include "lsst/meas/multifit/BaseEvaluator.h"

%shared_ptr(lsst::meas::multifit::Evaluator);

%include "lsst/meas/multifit/Evaluator.h"

%include "lsst/meas/multifit/Evaluation.h"

%{ 
#include "lsst/meas/multifit/SourceMeasurement.h"
%}


%include "lsst/meas/multifit/SourceMeasurement.h"

%define %Exposure(PIXTYPE)
    lsst::afw::image::Exposure<PIXTYPE, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>
%enddef

%extend lsst::meas::multifit::SourceMeasurement {
    %template(measure) lsst::meas::multifit::SourceMeasurement::measure< %Exposure(float) >;
    %template(measure) lsst::meas::multifit::SourceMeasurement::measure< %Exposure(double) >;
    %template(getModelImage) lsst::meas::multifit::SourceMeasurement::getModelImage< %Exposure(float) >;
    %template(getModelImage) lsst::meas::multifit::SourceMeasurement::getModelImage< %Exposure(double) >;
}
