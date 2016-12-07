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

%define modelfitLib_DOCSTRING
"
Basic routines to talk to lsst::meas::modelfit classes
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.meas.modelfit", docstring=modelfitLib_DOCSTRING) modelfitLib

%{
#include "lsst/afw/geom.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/geom/ellipses/PyPixelRegion.h"
#include "lsst/afw/table.h"
#include "lsst/afw/cameraGeom.h"
#include "lsst/afw/image.h"
#include "lsst/afw/math.h"
#include "lsst/afw/detection.h"
#include "lsst/shapelet.h"
#include "lsst/meas/base.h"
#include "lsst/meas/modelfit.h"
%}

%include "lsst/p_lsstSwig.i"
%initializeNumPy(meas_modelfit)
%{
#include "ndarray/swig.h"
#include "ndarray/converter/eigen.h"
%}
%include "ndarray.i"

%include "lsst/base.h"
%include "std_complex.i"

%lsst_exceptions();

%include "std_vector.i"
%import "lsst/afw/geom/geomLib.i"
%import "lsst/afw/geom/ellipses/ellipsesLib.i"
%import "lsst/afw/table/io/ioLib.i"
%import "lsst/afw/table/tableLib.i"
%import "lsst/afw/image/imageLib.i"
%import "lsst/afw/math/random.i"
%import "lsst/shapelet/shapeletLib.i"
%import "lsst/meas/base/baseLib.i"
%import "lsst/pex/config.h"

%template(EpochFootprintVector) std::vector<PTR(lsst::meas::modelfit::EpochFootprint)>;

%declareNumPyConverters(ndarray::Array<lsst::meas::modelfit::Scalar,1,0>);
%declareNumPyConverters(ndarray::Array<lsst::meas::modelfit::Scalar,1,1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::modelfit::Scalar,2,1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::modelfit::Scalar,2,2>);
%declareNumPyConverters(ndarray::Array<lsst::meas::modelfit::Scalar const,1,0>);
%declareNumPyConverters(ndarray::Array<lsst::meas::modelfit::Scalar const,1,1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::modelfit::Scalar const,2,1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::modelfit::Scalar const,2,2>);
%declareNumPyConverters(ndarray::Array<lsst::meas::modelfit::Scalar,2,-1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::modelfit::Scalar,2,-2>);
%declareNumPyConverters(ndarray::Array<lsst::meas::modelfit::Scalar const,2,-1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::modelfit::Scalar const,2,-2>);
%declareNumPyConverters(ndarray::Array<lsst::meas::modelfit::Pixel,1,0>);
%declareNumPyConverters(ndarray::Array<lsst::meas::modelfit::Pixel,1,1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::modelfit::Pixel,2,1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::modelfit::Pixel,2,2>);
%declareNumPyConverters(ndarray::Array<lsst::meas::modelfit::Pixel const,1,0>);
%declareNumPyConverters(ndarray::Array<lsst::meas::modelfit::Pixel const,1,1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::modelfit::Pixel const,2,1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::modelfit::Pixel const,2,2>);
%declareNumPyConverters(ndarray::Array<lsst::meas::modelfit::Pixel,2,-1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::modelfit::Pixel,2,-2>);
%declareNumPyConverters(ndarray::Array<lsst::meas::modelfit::Pixel const,2,-1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::modelfit::Pixel const,2,-2>);
%declareNumPyConverters(lsst::meas::modelfit::Vector);
%declareNumPyConverters(lsst::meas::modelfit::Matrix);
%declareNumPyConverters(Eigen::VectorXd);
%declareNumPyConverters(Eigen::MatrixXd);
%declareNumPyConverters(ndarray::Array<double,1,1>);
%declareNumPyConverters(ndarray::Array<double,2,2>);

%include "lsst/meas/modelfit/common.h"

%pythoncode %{
import numpy
Scalar = numpy.float64
Pixel = numpy.float32
%}

%shared_ptr(lsst::meas::modelfit::Prior);
%shared_ptr(lsst::meas::modelfit::MixturePrior);
%shared_ptr(lsst::meas::modelfit::SoftenedLinearPrior);
%shared_ptr(lsst::meas::modelfit::SemiEmpiricalPrior);

%shared_ptr(lsst::meas::modelfit::ModelFitTable);
%shared_ptr(lsst::meas::modelfit::ModelFitRecord);
%shared_ptr(lsst::meas::modelfit::Model);
%shared_ptr(lsst::meas::modelfit::MultiModel);
%shared_ptr(lsst::meas::modelfit::Interpreter);
%shared_ptr(lsst::meas::modelfit::Likelihood);
%shared_ptr(lsst::meas::modelfit::EpochFootprint);
%shared_ptr(lsst::meas::modelfit::UnitTransformedLikelihood);
%shared_ptr(lsst::meas::modelfit::Sampler);
%shared_ptr(lsst::meas::modelfit::SamplingObjective);
%shared_ptr(lsst::meas::modelfit::SamplingInterpreter);
%shared_ptr(lsst::meas::modelfit::DirectSamplingInterpreter);
%shared_ptr(lsst::meas::modelfit::MarginalSamplingInterpreter);
%shared_ptr(lsst::meas::modelfit::AdaptiveImportanceSampler);
%shared_ptr(lsst::meas::modelfit::MultiShapeletPsfLikelihood);

//----------- Mixtures --------------------------------------------------------------------------------------

%declareTablePersistable(Mixture, lsst::meas::modelfit::Mixture);
%ignore lsst::meas::modelfit::Mixture::begin;
%ignore lsst::meas::modelfit::Mixture::end;
%ignore lsst::meas::modelfit::Mixture::operator[];
%rename(__len__) lsst::meas::modelfit::Mixture::size;

%include "lsst/meas/modelfit/Mixture.h"

%ignore std::vector<lsst::meas::modelfit::MixtureComponent>::vector(size_type);
%ignore std::vector<lsst::meas::modelfit::MixtureComponent>::resize(size_type);
%template(MixtureComponentList) std::vector<lsst::meas::modelfit::MixtureComponent>;

%addStreamRepr(lsst::meas::modelfit::MixtureComponent);
%addStreamRepr(lsst::meas::modelfit::Mixture);

%extend lsst::meas::modelfit::Mixture {
    lsst::meas::modelfit::MixtureComponent & __getitem__(std::size_t i) {
        return (*($self))[i];
    }
    lsst::meas::modelfit::Scalar evaluate(
        lsst::meas::modelfit::MixtureComponent const & component,
        lsst::meas::modelfit::Vector const & x
    ) const {
        return $self->evaluate(component, x);
    }
    lsst::meas::modelfit::Scalar evaluate(
        lsst::meas::modelfit::Vector const & x
    ) const {
        return $self->evaluate(x);
    }

    %pythoncode %{
        def __iter__(self):
            for i in range(len(self)):
                yield self[i]
    %}
}

%pythoncode %{
    Mixture.UpdateRestriction = MixtureUpdateRestriction
    Mixture.Component = MixtureComponent
    Mixture.ComponentList = MixtureComponentList
%}

//----------- Miscellaneous ---------------------------------------------------------------------------------

%copyctor lsst::meas::modelfit::SoftenedLinearPriorControl;
%returnCopy(lsst::meas::modelfit::SoftenedLinearPrior::getControl)

%copyctor lsst::meas::modelfit::SemiEmpiricalPriorControl;
%returnCopy(lsst::meas::modelfit::SemiEmpiricalPrior::getControl)

%include "lsst/meas/modelfit/Model.h"
%include "lsst/meas/modelfit/MultiModel.h"
%include "lsst/meas/modelfit/Prior.h"
%include "lsst/meas/modelfit/MixturePrior.h"
%include "lsst/meas/modelfit/SoftenedLinearPrior.h"
%include "lsst/meas/modelfit/SemiEmpiricalPrior.h"
%include "lsst/meas/modelfit/Interpreter.h"
%include "lsst/meas/modelfit/Likelihood.h"
%include "lsst/meas/modelfit/UnitSystem.h"
%include "lsst/meas/modelfit/UnitTransformedLikelihood.h"
%include "lsst/meas/modelfit/Sampling.h"
%include "lsst/meas/modelfit/Sampler.h"
%include "lsst/meas/modelfit/DirectSamplingInterpreter.h"
%include "lsst/meas/modelfit/AdaptiveImportanceSampler.h"
%include "lsst/meas/modelfit/TruncatedGaussian.h"

// work around a bug in SWIG 3.0.2: mis-handling templated constructors
// once we have a fixed SWIG you may remove the hack from UnitSystem.h and uncomment the next few lines
// %extend lsst::meas::modelfit::UnitSystem {
//     %template(UnitSystem) UnitSystem<float>;
//     %template(UnitSystem) UnitSystem<double>;
// }

%pythoncode %{
import lsst.pex.config
SoftenedLinearPriorConfig = lsst.pex.config.makeConfigClass(SoftenedLinearPriorControl)
SemiEmpiricalPriorConfig = lsst.pex.config.makeConfigClass(SemiEmpiricalPriorControl)

SoftenedLinearPrior.Control = SoftenedLinearPriorControl
SemiEmpiricalPrior.Control = SemiEmpiricalPriorControl
SoftenedLinearPrior.ConfigClass = SoftenedLinearPriorConfig
SemiEmpiricalPrior.ConfigClass = SemiEmpiricalPriorConfig
%}

%castShared(lsst::meas::modelfit::SoftenedLinearPrior, lsst::meas::modelfit::Prior)
%castShared(lsst::meas::modelfit::SemiEmpiricalPrior, lsst::meas::modelfit::Prior)
%castShared(lsst::meas::modelfit::MixturePrior, lsst::meas::modelfit::Prior)
%castShared(lsst::meas::modelfit::MultiModel, lsst::meas::modelfit::Model)
%castShared(lsst::meas::modelfit::SamplingInterpreter, lsst::meas::modelfit::Interpreter)
%castShared(lsst::meas::modelfit::DirectSamplingInterpreter, lsst::meas::modelfit::Interpreter)
%castShared(lsst::meas::modelfit::MarginalSamplingInterpreter, lsst::meas::modelfit::Interpreter)

%ignore std::vector<lsst::afw::geom::ellipses::Ellipse>::vector(size_type);
%ignore std::vector<lsst::afw::geom::ellipses::Ellipse>::resize(size_type);
%template(EllipseVector) std::vector<lsst::afw::geom::ellipses::Ellipse>;
%template(NameVector) std::vector<std::string>;
%template(BasisVector) std::vector<PTR(lsst::shapelet::MultiShapeletBasis)>;
%template(ModelVector) std::vector<PTR(lsst::meas::modelfit::Model)>;

%pythoncode %{
Model.EllipseVector = EllipseVector
Model.BasisVector = BasisVector
Model.NameVector = NameVector
%}

%include "std_map.i"
%template(ImportanceSamplerControlMap) std::map<int,lsst::meas::modelfit::ImportanceSamplerControl>;

%extend lsst::meas::modelfit::AdaptiveImportanceSampler {
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

UnitTransformedLikelihoodConfig = lsst.pex.config.makeConfigClass(UnitTransformedLikelihoodControl)
UnitTransformedLikelihood.ConfigClass = UnitTransformedLikelihoodConfig
%}

//----------- ModelFitRecord/Table/Catalog ------------------------------------------------------------------

%include "lsst/meas/modelfit/ModelFitRecord.h"

%addCastMethod(lsst::meas::modelfit::ModelFitTable, lsst::afw::table::BaseTable)
%addCastMethod(lsst::meas::modelfit::ModelFitRecord, lsst::afw::table::BaseRecord)


%template(ModelFitColumnView) lsst::afw::table::ColumnViewT<lsst::meas::modelfit::ModelFitRecord>;

%include "lsst/afw/table/SortedCatalog.i"

namespace lsst { namespace afw { namespace table {

using meas::modelfit::ModelFitRecord;
using meas::modelfit::ModelFitTable;

%declareSortedCatalog(SortedCatalogT, ModelFit)

}}} // namespace lsst::afw::table

namespace lsst { namespace meas { namespace modelfit {

typedef lsst::afw::table::SortedCatalogT<ModelFitRecord> ModelFitCatalog;

}}} // namespace lsst::meas::modelfit

//----------- More Miscellaneous ----------------------------------------------------------------------------

%include "lsst/meas/modelfit/integrals.h"
%include "lsst/meas/modelfit/optimizer.i"
%include "lsst/meas/modelfit/MarginalSamplingInterpreter.h"
%include "lsst/meas/modelfit/psf.i"
%include "lsst/meas/modelfit/CModel.i"
