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

#include "ndarray/eigen.h"

#include "lsst/meas/multifit/MarginalSampling.h"
#include "lsst/meas/multifit/ModelFitRecord.h"

namespace lsst { namespace meas { namespace multifit {

namespace {

class MarginalSamplingObjective : public SamplingObjective {
public:

    virtual Scalar operator()(
        ndarray::Array<Scalar const,1,1> const & parameters,
        afw::table::BaseRecord & sample
    ) const {
        _likelihood->computeModelMatrix(_modelMatrix, parameters);
        Matrix hessian = Matrix::Zero(_likelihood->getAmplitudeDim(), _likelihood->getAmplitudeDim());
        hessian.selfadjointView<Eigen::Lower>().rankUpdate(_modelMatrix.asEigen().adjoint().cast<Scalar>());
        Vector gradient
            = -(_modelMatrix.asEigen().adjoint() * _likelihood->getData().asEigen()).cast<Scalar>();
        ArrayKey nestedKey = static_cast<MarginalSamplingInterpreter &>(*getInterpreter()).getNestedKey();
        ndarray::Array<Scalar,1,1> nested = sample[nestedKey];
        int const n = _likelihood->getAmplitudeDim();
        for (int i = 0, k = n; i < n; ++i) {
            nested[i] = gradient[i];
            for (int j = 0; j <= i; ++j, ++k) {
                nested[k] = hessian(i, j);
            }
        }
        return getInterpreter()->getPrior()->marginalize(gradient, hessian, parameters);
    }

    MarginalSamplingObjective(
        PTR(SamplingInterpreter) interpreter,
        PTR(Likelihood) likelihood
    ) : SamplingObjective(interpreter, likelihood)
    {
        if (!getInterpreter()->getPrior()) {
            throw LSST_EXCEPT(
                pex::exceptions::LogicErrorException,
                "Cannot create MarginalSamplingObjective without a Prior"
            );
        }
    }

};

} // anonymous

MarginalSamplingInterpreter::MarginalSamplingInterpreter(
    afw::table::Schema & sampleSchema,
    PTR(Model) model,
    PTR(Prior) prior
) :
    SamplingInterpreter(
        sampleSchema,
        model->getNonlinearNames(),
        model,
        prior
    )
{
    _parameterKey = sampleSchema.addField(
        afw::table::Field< afw::table::Array<Scalar> >(
            "parameters", "nonlinear parameters at this sample point (amplitudes are nested)",
            model->getNonlinearDim()
        ),
        true // doReplace
    );
    int n = model->getAmplitudeDim();
    _nestedKey = sampleSchema.addField(
        afw::table::Field< afw::table::Array<Scalar> >(
            "nested", "an opaque representation of the nested amplitude likelihood",
            n + n*(n+1)/2
        ),
        true
    );
    _nonlinearKey = _parameterKey;
}

ndarray::Array<Scalar,1,1> MarginalSamplingInterpreter::computeAmplitudeQuantiles(
    ModelFitRecord const & record,
    ndarray::Array<Scalar const,1,1> const & fractions,
    int index
) const {
    throw LSST_EXCEPT(
        pex::exceptions::LogicErrorException,
        "Amplitude quantiles not implemented for marginalized samples"
    );
}

ndarray::Array<Scalar,1,1> MarginalSamplingInterpreter::computeAmplitudeMean(
    ModelFitRecord const & record
) const {
    throw LSST_EXCEPT(
        pex::exceptions::LogicErrorException,
        "Amplitude means not implemented for marginalized samples"
    );
}

ndarray::Array<Scalar,2,2> MarginalSamplingInterpreter::computeAmplitudeCovariance(
    ModelFitRecord const & record,
    ndarray::Array<Scalar const,1,1> const & mean
) const {
    throw LSST_EXCEPT(
        pex::exceptions::LogicErrorException,
        "Amplitude covariances not implemented for marginalized samples"
    );
}

PTR(SamplingObjective) MarginalSamplingInterpreter::makeObjective(
    PTR(SamplingInterpreter) self,
    PTR(Likelihood) likelihood
) const {
    return boost::make_shared<MarginalSamplingObjective>(self, likelihood);
}

PTR(Interpreter) MarginalSamplingInterpreter::_clone() const {
    return boost::make_shared<MarginalSamplingInterpreter>(*this);
}

void MarginalSamplingInterpreter::_packParameters(
    ndarray::Array<Scalar const,1,1> const & nonlinear,
    ndarray::Array<Scalar const,1,1> const & amplitudes,
    ndarray::Array<Scalar,1,1> const & parameters
) const {
    parameters.deep() = nonlinear;
}

void MarginalSamplingInterpreter::_unpackNonlinear(
    ndarray::Array<Scalar const,1,1> const & parameters,
    ndarray::Array<Scalar,1,1> const & nonlinear
) const {
    nonlinear.deep() = parameters;
}

ModelFitCatalog MarginalSamplingInterpreter::unnset(
    ModelFitCatlog const & inCat,
    afw::math::Random & rng,
    double multiplier
) {
    afw::table::SchemaMapper mapper(inCat.getTable()->getSampleTable()->getSchema());
    
}

}}} // namespace lsst::meas::multifit
