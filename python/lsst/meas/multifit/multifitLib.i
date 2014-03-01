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
#include "lsst/shapelet.h"
#include "lsst/meas/algorithms.h"
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

%include "std_vector.i"
%include "ndarray.i"
%import "lsst/afw/geom/geomLib.i"
%import "lsst/afw/geom/ellipses/ellipsesLib.i"
%import "lsst/afw/table/io/ioLib.i"
%import "lsst/afw/table/tableLib.i"
%import "lsst/afw/image/imageLib.i"
%import "lsst/afw/math/random.i"
%import "lsst/meas/algorithms/algorithmsLib.i"
%import "lsst/shapelet/shapeletLib.i"
%import "lsst/pex/config.h"

%template(EpochFootprintVector) std::vector<PTR(lsst::meas::multifit::EpochFootprint)>;

%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Scalar,1,0>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Scalar,1,1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Scalar,2,1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Scalar,2,2>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Scalar const,1,0>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Scalar const,1,1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Scalar const,2,1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Scalar const,2,2>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Scalar,2,-1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Scalar,2,-2>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Scalar const,2,-1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Scalar const,2,-2>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Pixel,1,0>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Pixel,1,1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Pixel,2,1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Pixel,2,2>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Pixel const,1,0>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Pixel const,1,1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Pixel const,2,1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Pixel const,2,2>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Pixel,2,-1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Pixel,2,-2>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Pixel const,2,-1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Pixel const,2,-2>);
%declareNumPyConverters(lsst::meas::multifit::Vector);
%declareNumPyConverters(lsst::meas::multifit::Matrix);
%declareNumPyConverters(Eigen::VectorXd);
%declareNumPyConverters(Eigen::MatrixXd);
%declareNumPyConverters(ndarray::Array<double,1,1>);
%declareNumPyConverters(ndarray::Array<double,2,2>);

%include "lsst/meas/multifit/constants.h"

%pythoncode %{
import numpy
Scalar = numpy.float64
Pixel = numpy.float32
%}

%declareTablePersistable(Prior, lsst::meas::multifit::Prior);
%declareTablePersistable(MixturePrior, lsst::meas::multifit::MixturePrior);

%shared_ptr(lsst::meas::multifit::ModelFitTable);
%shared_ptr(lsst::meas::multifit::ModelFitRecord);
%shared_ptr(lsst::meas::multifit::Model);
%shared_ptr(lsst::meas::multifit::MultiModel);
%shared_ptr(lsst::meas::multifit::Interpreter);
%shared_ptr(lsst::meas::multifit::Likelihood);
%shared_ptr(lsst::meas::multifit::EpochFootprint);
%shared_ptr(lsst::meas::multifit::ProjectedLikelihood);
%shared_ptr(lsst::meas::multifit::Sampler);
%shared_ptr(lsst::meas::multifit::SamplingObjective);
%shared_ptr(lsst::meas::multifit::SamplingInterpreter);
%shared_ptr(lsst::meas::multifit::DirectSamplingInterpreter);
%shared_ptr(lsst::meas::multifit::MarginalSamplingInterpreter);
%shared_ptr(lsst::meas::multifit::AdaptiveImportanceSampler);
%shared_ptr(lsst::meas::multifit::MultiShapeletPsfLikelihood);

//----------- Mixtures --------------------------------------------------------------------------------------

%declareTablePersistable(Mixture, lsst::meas::multifit::Mixture);
%ignore lsst::meas::multifit::Mixture::begin;
%ignore lsst::meas::multifit::Mixture::end;
%ignore lsst::meas::multifit::Mixture::operator[];
%rename(__len__) lsst::meas::multifit::Mixture::size;

%include "lsst/meas/multifit/Mixture.h"

%ignore std::vector<lsst::meas::multifit::MixtureComponent>::vector(size_type);
%ignore std::vector<lsst::meas::multifit::MixtureComponent>::resize(size_type);
%template(MixtureComponentList) std::vector<lsst::meas::multifit::MixtureComponent>;

%addStreamRepr(lsst::meas::multifit::MixtureComponent);
%addStreamRepr(lsst::meas::multifit::Mixture);

%extend lsst::meas::multifit::Mixture {
    lsst::meas::multifit::MixtureComponent & __getitem__(std::size_t i) {
        return (*($self))[i];
    }
    lsst::meas::multifit::Scalar evaluate(
        lsst::meas::multifit::MixtureComponent const & component,
        lsst::meas::multifit::Vector const & x
    ) const {
        return $self->evaluate(component, x);
    }
    lsst::meas::multifit::Scalar evaluate(
        lsst::meas::multifit::Vector const & x
    ) const {
        return $self->evaluate(x);
    }

    %pythoncode %{
        def __iter__(self):
            for i in xrange(len(self)):
                yield self[i]
    %}
}

%pythoncode %{
    Mixture.UpdateRestriction = MixtureUpdateRestriction
    Mixture.Component = MixtureComponent
    Mixture.ComponentList = MixtureComponentList
%}

//----------- Miscellaneous ---------------------------------------------------------------------------------

%include "lsst/meas/multifit/models.h"
%include "lsst/meas/multifit/priors.h"
%include "lsst/meas/multifit/Interpreter.h"
%include "lsst/meas/multifit/Likelihood.h"
%include "lsst/meas/multifit/UnitSystem.h"
%include "lsst/meas/multifit/ProjectedLikelihood.h"
%include "lsst/meas/multifit/Sampling.h"
%include "lsst/meas/multifit/DirectSampling.h"
%include "lsst/meas/multifit/MarginalSampling.h"
%include "lsst/meas/multifit/AdaptiveImportanceSampler.h"
%include "lsst/meas/multifit/psf.h"
%include "lsst/meas/multifit/TruncatedGaussian.h"

%extend lsst::meas::multifit::UnitSystem {
    %template(UnitSystem) UnitSystem<float>;
    %template(UnitSystem) UnitSystem<double>;
}

%castShared(lsst::meas::multifit::MixturePrior, lsst::meas::multifit::Prior)
%castShared(lsst::meas::multifit::MultiModel, lsst::meas::multifit::Model)
%castShared(lsst::meas::multifit::SamplingInterpreter, lsst::meas::multifit::Interpreter)
%castShared(lsst::meas::multifit::DirectSamplingInterpreter, lsst::meas::multifit::Interpreter)
%castShared(lsst::meas::multifit::MarginalSamplingInterpreter, lsst::meas::multifit::Interpreter)

%ignore std::vector<lsst::afw::geom::ellipses::Ellipse>::vector(size_type);
%ignore std::vector<lsst::afw::geom::ellipses::Ellipse>::resize(size_type);
%template(EllipseVector) std::vector<lsst::afw::geom::ellipses::Ellipse>;
%template(NameVector) std::vector<std::string>;
%template(BasisVector) std::vector<PTR(lsst::shapelet::MultiShapeletBasis)>;
%template(ModelVector) std::vector<PTR(lsst::meas::multifit::Model)>;

%pythoncode %{
Model.EllipseVector = EllipseVector
Model.BasisVector = BasisVector
Model.NameVector = NameVector
%}

%include "std_map.i"
%template(ImportanceSamplerControlMap) std::map<int,lsst::meas::multifit::ImportanceSamplerControl>;

%extend lsst::meas::multifit::AdaptiveImportanceSampler {
%pythoncode %{
@property
def iterations(self):
    d = {}
    for k, v in self.getIterations().items():
        for i, s in enumerate(v):
            d[k,i] = s
    return d
%}
}

%pythoncode %{
import lsst.pex.config
import numpy

ProjectedLikelihoodConfig = lsst.pex.config.makeConfigClass(ProjectedLikelihoodControl)
ProjectedLikelihood.ConfigClass = ProjectedLikelihoodConfig
%}

//----------- ModelFitRecord/Table/Catalog ------------------------------------------------------------------

%include "lsst/meas/multifit/ModelFitRecord.h"

%addCastMethod(lsst::meas::multifit::ModelFitTable, lsst::afw::table::BaseTable)
%addCastMethod(lsst::meas::multifit::ModelFitRecord, lsst::afw::table::BaseRecord)


%template(ModelFitColumnView) lsst::afw::table::ColumnViewT<lsst::meas::multifit::ModelFitRecord>;

%include "lsst/afw/table/SortedCatalog.i"

namespace lsst { namespace afw { namespace table {

using meas::multifit::ModelFitRecord;
using meas::multifit::ModelFitTable;

%declareSortedCatalog(SortedCatalogT, ModelFit)

}}} // namespace lsst::afw::table

namespace lsst { namespace meas { namespace multifit {

typedef lsst::afw::table::SortedCatalogT<ModelFitRecord> ModelFitCatalog;

}}} // namespace lsst::meas::multifit

//----------- More Miscellaneous ----------------------------------------------------------------------------

%include "lsst/meas/multifit/integrals.h"
%include "lsst/meas/multifit/optimizer.i"
%include "lsst/meas/multifit/CModel.i"
