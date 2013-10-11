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

#include "lsst/afw/table/BaseRecord.h"
#include "lsst/meas/multifit/Sampler.h"

namespace lsst { namespace meas { namespace multifit {

namespace {

class DirectSamplerObjective : public SamplerObjective {
public:

    virtual Scalar operator()(
        ndarray::Array<Scalar const,1,1> const & parameters,
        afw::table::BaseRecord & sample
    ) const {
        int np = _likelihood->getNonlinearDim();
        ndarray::Array<Scalar const,1,1> nonlinear
            = parameters[ndarray::view(0, np)];
        ndarray::Array<Scalar const,1,1> amplitudes
            = parameters[ndarray::view(np, np + _likelihood->getAmplitudeDim())];
        _likelihood->computeModelMatrix(_modelMatrix, nonlinear);
        Scalar chiSq = (
            _likelihood->getData().asEigen() - _modelMatrix.asEigen() * amplitudes.asEigen().cast<Pixel>()
        ).squaredNorm();
        if (_prior) {
            chiSq -= std::log(_prior->evaluate(nonlinear, amplitudes));
        }
        return chiSq;
    }

    DirectSamplerObjective(
        ArrayKey const & parameterKey,
        PTR(Likelihood) likelihood,
        PTR(Prior) prior
    ) :
        SamplerObjective(parameterKey, likelihood, prior)
    {}

private:
    ndarray::Array<Pixel,1,1> _residuals;
};

class MarginalizingSamplerObjective : public SamplerObjective {
public:

    virtual Scalar operator()(
        ndarray::Array<Scalar const,1,1> const & parameters,
        afw::table::BaseRecord & sample
    ) const {
        _likelihood->computeModelMatrix(_modelMatrix, parameters);
        Matrix fisher = Matrix::Zero(_likelihood->getAmplitudeDim(), _likelihood->getAmplitudeDim());
        fisher.selfadjointView<Eigen::Lower>().rankUpdate(_modelMatrix.asEigen().adjoint().cast<Scalar>());
        Vector gradient
            = -(_modelMatrix.asEigen().adjoint() * _likelihood->getData().asEigen()).cast<Scalar>();
        ndarray::Array<Scalar,1,1> nested = sample[_nestedKey];
        int const n = _likelihood->getAmplitudeDim();
        for (int i = 0, k = n; i < n; ++i) {
            nested[i] = gradient[i];
            for (int j = 0; j <= i; ++j, ++k) {
                nested[k] = fisher(i, j);
            }
        }
        return _prior->marginalize(gradient, fisher, parameters);
    }

    MarginalizingSamplerObjective(
        ArrayKey const & parameterKey,
        ArrayKey const & nestedKey,
        PTR(Likelihood) likelihood,
        PTR(Prior) prior
    ) : SamplerObjective(parameterKey, likelihood, prior),
        _nestedKey(nestedKey)
    {
        if (!_prior) {
            throw LSST_EXCEPT(
                pex::exceptions::LogicErrorException,
                "Cannot create marginalizing SamplerObjective without a Prior"
            );
        }
    }

private:
    ArrayKey _nestedKey;
};

} // anonymous

SamplerObjective::SamplerObjective(
    ArrayKey const & parameterKey,
    PTR(Likelihood) likelihood,
    PTR(Prior) prior
) :
    _parameterKey(parameterKey),
    _likelihood(likelihood),
    _prior(prior),
    _modelMatrix(ndarray::allocate(likelihood->getDataDim(), likelihood->getAmplitudeDim()))
{}

PTR(SamplerObjective) SamplerObjectiveFactory::operator()(PTR(Likelihood) likelihood) const {
    if (_doMarginalizeAmplitudes) {
        return boost::make_shared<MarginalizingSamplerObjective>(
            _parameterKey, _nestedKey, likelihood, _prior
        );
    } else {
        return boost::make_shared<DirectSamplerObjective>(_parameterKey, likelihood, _prior);
    }
}

void SamplerObjectiveFactory::mapParameters(
    ndarray::Array<Scalar const,1,1> const & nonlinear,
    ndarray::Array<Scalar const,1,1> const & amplitudes,
    ndarray::Array<Scalar,1,1> const & parameters
) const {
    LSST_ASSERT_EQUAL(
        parameters.getSize<0>(), getParameterDim(),
        "Size of parameter array (%d) does not match expected size (%d)",
        pex::exceptions::LengthErrorException
    );
    if (_doMarginalizeAmplitudes) {
        LSST_ASSERT_EQUAL(
            nonlinear.getSize<0>(), getParameterDim(),
            "Size of nonlinear array (%d) does not match expected size (%d)",
            pex::exceptions::LengthErrorException
        );
        parameters.deep() = nonlinear;
    } else {
        int n = nonlinear.getSize<0>();
        LSST_ASSERT_EQUAL(
            n + amplitudes.getSize<0>(), getParameterDim(),
            "Combined size of nonlinear and amplitude arrays (%d) does not match expected size (%d)",
            pex::exceptions::LengthErrorException
        );
        parameters[ndarray::view(0, n)] = nonlinear;
        parameters[ndarray::view(n, n+amplitudes.getSize<0>())] = amplitudes;
    }
}

SamplerObjectiveFactory::SamplerObjectiveFactory(
    afw::table::Schema & sampleSchema,
    PTR(Model) model,
    PTR(Prior) prior,
    bool doMarginalizeAmplitudes
) : _doMarginalizeAmplitudes(doMarginalizeAmplitudes), _prior(prior) {
    if (doMarginalizeAmplitudes) {
        _parameterKey = sampleSchema.addField< afw::table::Array<Scalar> >(
            "parameters", "nonlinear parameters at this sample point (amplitudes are nested)",
            model->getNonlinearDim()
        );
        int n = model->getAmplitudeDim();
        // TODO: need to provide public interface for interpreting the nested distribution
        _nestedKey = sampleSchema.addField< afw::table::Array<Scalar> >(
            "nested", "an opaque representation of the nested amplitude distribution",
            n + n*(n+1)/2
        );
    } else {
        _parameterKey = sampleSchema.addField< afw::table::Array<Scalar> >(
            "parameters", "nonlinear and amplitude parameters at this sample point",
            model->getNonlinearDim() + model->getAmplitudeDim()
        );
    }
}

}}} // namespace lsst::meas::multifit
