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
#include "lsst/pex/logging.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/geom/ellipses/PyPixelRegion.h"
#include "lsst/afw/table.h"
#include "lsst/afw/cameraGeom.h"
#include "lsst/afw/image.h"
#include "lsst/meas/multifit.h"
#include "ndarray/eigen.h"
#include "Eigen/Core"

// namespace-ish hack required by NumPy C-API; see NumPy docs for more info
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
%import "lsst/afw/table/io/ioLib.i"
%import "lsst/afw/table/tableLib.i"
%import "lsst/afw/image/imageLib.i"
%import "lsst/pex/config.h"

namespace lsst { namespace shapelet {
class MultiShapeletBasis;
}}

%declareNumPyConverters(lsst::meas::multifit::Vector);
%declareNumPyConverters(lsst::meas::multifit::Matrix);
%declareNumPyConverters(Eigen::VectorXd);
%declareNumPyConverters(Eigen::MatrixXd);
%declareNumPyConverters(ndarray::Array<double,1,1>);
%declareNumPyConverters(ndarray::Array<double,2,2>);

%declareTablePersistable(SampleSet, lsst::meas::multifit::SampleSet);

%shared_ptr(lsst::meas::multifit::Prior);
%shared_ptr(lsst::meas::multifit::SingleComponentPrior);
%shared_ptr(lsst::meas::multifit::ExpectationFunctor);
%shared_ptr(lsst::meas::multifit::Objective);
%shared_ptr(lsst::meas::multifit::SingleEpochObjective);
%shared_ptr(lsst::meas::multifit::BaseSampler);
%shared_ptr(lsst::meas::multifit::NaiveGridSampler);

%include "lsst/meas/multifit/constants.h"
%include "lsst/meas/multifit/LogGaussian.h"
%include "lsst/meas/multifit/priors.h"
%include "lsst/meas/multifit/Objective.h"
%include "lsst/meas/multifit/KernelDensityEstimator.h"
%include "lsst/meas/multifit/SampleSet.h"
%include "lsst/meas/multifit/ExpectationFunctor.h"
%include "lsst/meas/multifit/BaseSampler.h"
%include "lsst/meas/multifit/NaiveGridSampler.h"

%pythoncode %{
import lsst.pex.config
import numpy
SingleEpochObjectiveConfig = lsst.pex.config.makeConfigClass(SingleEpochObjectiveControl)
SingleEpochObjective.ConfigClass = SingleEpochObjectiveConfig
%}

%shared_ptr(lsst::meas::multifit::ModelFitTable);
%shared_ptr(lsst::meas::multifit::ModelFitRecord);

%include "lsst/meas/multifit/tables.h"

%addCastMethod(lsst::meas::multifit::ModelFitTable, lsst::afw::table::BaseTable)
%addCastMethod(lsst::meas::multifit::ModelFitRecord, lsst::afw::table::BaseRecord)


%template(ModelFitColumnView) lsst::afw::table::ColumnViewT<lsst::meas::multifit::ModelFitRecord>;

%include "lsst/afw/table/SortedCatalog.i"

namespace lsst { namespace afw { namespace table {

using meas::multifit::ModelFitRecord;
using meas::multifit::ModelFitTable;

%declareSortedCatalog(SortedCatalogT, ModelFit)

}}} // namespace lsst::afw::table
