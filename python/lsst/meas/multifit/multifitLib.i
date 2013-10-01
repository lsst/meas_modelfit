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
%import "lsst/shapelet/shapeletLib.i"
%import "lsst/pex/config.h"

%template(VectorEpochFootprint) std::vector<PTR(lsst::meas::multifit::EpochFootprint)>;

%declareNumPyConverters(lsst::meas::multifit::samples::Vector);
%declareNumPyConverters(lsst::meas::multifit::samples::Matrix);
%declareNumPyConverters(Eigen::VectorXd);
%declareNumPyConverters(Eigen::MatrixXd);
%declareNumPyConverters(ndarray::Array<double,1,1>);
%declareNumPyConverters(ndarray::Array<double,2,2>);

%include "lsst/meas/multifit/constants.h"

%declareTablePersistable(Prior, lsst::meas::multifit::Prior);
%declareTablePersistable(FlatPrior, lsst::meas::multifit::FlatPrior);
%declareTablePersistable(MixturePrior, lsst::meas::multifit::MixturePrior);

%shared_ptr(lsst::meas::multifit::ParameterConverter);
%shared_ptr(lsst::meas::multifit::ParameterDefinition);
%shared_ptr(lsst::meas::multifit::Model);
%shared_ptr(lsst::meas::multifit::MultiModel);
%shared_ptr(lsst::meas::multifit::Likelihood);
%shared_ptr(lsst::meas::multifit::SingleEpochLikelihood);
%shared_ptr(lsst::meas::multifit::EpochFootprint);
%shared_ptr(lsst::meas::multifit::MultiEpochLikelihood);
%shared_ptr(lsst::meas::multifit::Sampler);
%shared_ptr(lsst::meas::multifit::AdaptiveImportanceSampler);

//----------- Mixtures --------------------------------------------------------------------------------------

namespace lsst { namespace meas { namespace multifit {

template <int N> class Mixture;
template <int N> class MixtureUpdateRestriction;
template <int N> class MixtureComponent;

}}} // namespace lsst::meas::multifit

%declareTablePersistable(MixtureBase, lsst::meas::multifit::MixtureBase);

%include "lsst/meas/multifit/MixtureBase.h"
%include "lsst/meas/multifit/Mixture.h"

%pythoncode %{
    Mixture = dict()
    MixtureComponent = dict()
%}

%addStreamRepr(lsst::meas::multifit::MixtureComponent);
%addStreamRepr(lsst::meas::multifit::MixtureBase);

%define %instantiateMixture(N)
%declareTablePersistable(Mixture ## N, lsst::meas::multifit::Mixture<N>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Scalar,1,0>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Scalar,1,1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Scalar,2,1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Scalar const,1,0>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Scalar const,1,1>);
%declareNumPyConverters(ndarray::Array<lsst::meas::multifit::Scalar const,2,1>);
%declareNumPyConverters(lsst::meas::multifit::MixtureComponent<N>::Vector);
%declareNumPyConverters(lsst::meas::multifit::MixtureComponent<N>::Matrix);
%template(MixtureComponent ## N) lsst::meas::multifit::MixtureComponent<N>;
%template(MixtureUpdateRestriction ## N) lsst::meas::multifit::MixtureUpdateRestriction<N>;
%template(MixtureComponent ## N ## List) std::vector<
    lsst::meas::multifit::Mixture<N>::Component,
    Eigen::aligned_allocator<lsst::meas::multifit::Mixture<N>::Component>
    >;
%ignore lsst::meas::multifit::Mixture<N>::begin;
%ignore lsst::meas::multifit::Mixture<N>::end;
%ignore lsst::meas::multifit::Mixture<N>::operator[];
%rename(__len__) lsst::meas::multifit::Mixture<N>::size;
%extend lsst::meas::multifit::Mixture<N> {
    lsst::meas::multifit::MixtureComponent<N> & __getitem__(std::size_t i) {
        return (*($self))[i];
    }
    lsst::meas::multifit::Scalar evaluate(
        lsst::meas::multifit::MixtureComponent<N> const & component,
        lsst::meas::multifit::MixtureComponent<N>::Vector const & x
    ) const {
        return $self->evaluate(component, x);
    }
    lsst::meas::multifit::Scalar evaluate(
        lsst::meas::multifit::MixtureComponent<N>::Vector const & x
    ) const {
        return $self->evaluate(x);
    }
    static PTR(lsst::meas::multifit::Mixture<N>) cast(PTR(lsst::meas::multifit::MixtureBase) const & p) {
        return boost::dynamic_pointer_cast< lsst::meas::multifit::Mixture<N> >(p);
    }
    %pythoncode %{
        UpdateRestriction = MixtureUpdateRestriction##N
        Component = MixtureComponent##N
        ComponentList = MixtureComponent##N##List
        def __iter__(self):
            for i in xrange(len(self)):
                yield self[i]
    %}
}
%template(Mixture ## N) lsst::meas::multifit::Mixture<N>;
%pythoncode %{
MixtureComponent[N] = MixtureComponent ## N
Mixture[N] = Mixture ## N
%}
%enddef

%instantiateMixture(1)
%instantiateMixture(2)
%instantiateMixture(3)

//----------- Miscellaneous ---------------------------------------------------------------------------------

%include "lsst/meas/multifit/LogGaussian.h"
%include "lsst/meas/multifit/parameters.h"
%include "lsst/meas/multifit/models.h"
%include "lsst/meas/multifit/priors.h"
%include "lsst/meas/multifit/Likelihood.h"
%include "lsst/meas/multifit/SingleEpochLikelihood.h"
%include "lsst/meas/multifit/MultiEpochLikelihood.h"
%include "lsst/meas/multifit/Sampler.h"
%include "lsst/meas/multifit/AdaptiveImportanceSampler.h"

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

SingleEpochLikelihoodConfig = lsst.pex.config.makeConfigClass(SingleEpochLikelihoodControl)
SingleEpochLikelihood.ConfigClass = SingleEpochLikelihoodConfig

MultiEpochLikelihoodConfig = lsst.pex.config.makeConfigClass(MultiEpochLikelihoodControl)
MultiEpochLikelihood.ConfigClass = MultiEpochLikelihoodConfig
%}

//----------- ModelFitRecord/Table/Catalog ------------------------------------------------------------------

%shared_ptr(lsst::meas::multifit::ModelFitTable);
%shared_ptr(lsst::meas::multifit::ModelFitRecord);

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


%include "lsst/meas/multifit/integrals.h"

%include "lsst/meas/multifit/optimizer.i"
