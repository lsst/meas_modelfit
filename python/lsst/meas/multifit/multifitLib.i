// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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

%{
#include "lsst/afw/geom.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/geom/ellipses/PyPixelRegion.h"
#include "lsst/meas/multifit.h"
#include "ndarray/eigen.h"
#include "Eigen/Core"
#define PY_ARRAY_UNIQUE_SYMBOL LSST_MEAS_MULTIFIT_NUMPY_ARRAY_API
#include "numpy/arrayobject.h"
#include "ndarray/swig.h"
#include "ndarray/swig/eigen.h"
%}

%init %{
    import_array();
%}

%include "lsst/p_lsstSwig.i"
%include "lsst/base.h"
%include "std_complex.i"

%lsst_exceptions();

%include "ndarray.i"
%import "lsst/afw/geom/geomLib.i"
%import "lsst/afw/geom/ellipses/ellipsesLib.i"
%import "lsst/pex/config.h"

namespace lsst { namespace shapelet {
class MultiShapeletBasis;
}}

%shared_ptr(lsst::meas::multifit::Prior);
%shared_ptr(lsst::meas::multifit::Objective);
%shared_ptr(lsst::meas::multifit::SingleEpochObjective);
%shared_ptr(lsst::meas::multifit::BaseSampler);

%include "lsst/meas/multifit/constants.h"
%include "lsst/meas/multifit/LogGaussian.h"
%include "lsst/meas/multifit/priors.h"
%include "lsst/meas/multifit/Objective.h"
%include "lsst/meas/multifit/BaseSampler.h"
