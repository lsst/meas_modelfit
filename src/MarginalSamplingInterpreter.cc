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

#include "boost/math/special_functions/round.hpp"

#include "ndarray/eigen.h"

#include "lsst/meas/multifit/MarginalSamplingInterpreter.h"
#include "lsst/meas/multifit/DirectSamplingInterpreter.h"
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

void MarginalSamplingInterpreter::unpackNested(
    ndarray::Array<Scalar const,1,1> const & nested, Vector & gradient, Matrix & hessian
) const {
    int const n = getAmplitudeDim();
    for (int i = 0, k = n; i < n; ++i) {
        gradient[i] = nested[i];
        for (int j = 0; j <= i; ++j, ++k) {
            hessian(i, j) = hessian(j, i) = nested[k];
        }
    }
}

void MarginalSamplingInterpreter::unpackNested(
    ndarray::Array<Scalar const,1,1> const & nested,
    ndarray::Array<Scalar,1,1> const & gradient,
    ndarray::Array<Scalar,2,2> const & hessian
) const {
    int const n = getAmplitudeDim();
    for (int i = 0, k = n; i < n; ++i) {
        gradient[i] = nested[i];
        for (int j = 0; j <= i; ++j, ++k) {
            hessian[i][j] = hessian[j][i] = nested[k];
        }
    }
}

ndarray::Array<Scalar,1,1> MarginalSamplingInterpreter::computeAmplitudeQuantiles(
    ModelFitRecord const & record,
    ndarray::Array<Scalar const,1,1> const & fractions,
    int index
) const {
    ndarray::Array<Scalar,1,1> result = ndarray::allocate(getAmplitudeDim());
    result.deep() = std::numeric_limits<Scalar>::quiet_NaN();
    return result;
}

ndarray::Array<Scalar,1,1> MarginalSamplingInterpreter::computeAmplitudeMean(
    ModelFitRecord const & record
) const {
    ndarray::Array<Scalar,1,1> result = ndarray::allocate(getAmplitudeDim());
    result.deep() = std::numeric_limits<Scalar>::quiet_NaN();
    return result;
}

ndarray::Array<Scalar,2,2> MarginalSamplingInterpreter::computeAmplitudeCovariance(
    ModelFitRecord const & record,
    ndarray::Array<Scalar const,1,1> const & mean
) const {
    ndarray::Array<Scalar,2,2> result = ndarray::allocate(getAmplitudeDim(), getAmplitudeDim());
    result.deep() = std::numeric_limits<Scalar>::quiet_NaN();
    return result;
}

PTR(SamplingObjective) MarginalSamplingInterpreter::makeObjective(
    PTR(SamplingInterpreter) self,
    PTR(Likelihood) likelihood
) const {
    return boost::make_shared<MarginalSamplingObjective>(self, likelihood);
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

namespace {

// This functor is applied to every field in a direct sampling schema; it
// finds the corresponding field in a marginal sampling schema (the input
// schema of the mapper), and thus creates a mapper that reproduces the
// direct schema while transferring appropriate fields from the marginal
// schcema.
struct MapSamplingFields {

    template <typename T>
    void operator()(afw::table::SchemaItem<T> const & directItem) const {
        afw::table::Key<T> directKey;
        if (directItem.field.getName() == "parameters") {
            // this is the only field we expect to be different
            directKey = mapper->addOutputField(directItem.field);
        } else {
            afw::table::Key<T> marginalKey = mapper->getInputSchema()[directItem.field.getName()];
            directKey = mapper->addMapping(marginalKey);
        }
        if (directKey != directItem.key) {
            throw LSST_EXCEPT(
                pex::exceptions::LogicErrorException,
                (boost::format("Unexpected mismatch between direct and marginal schemas at field '%s'")
                 % directItem.field.getName()).str()
            );
        }
    };

    afw::table::SchemaMapper * mapper;
};

} // anonymous


UnnestMarginalSamples::UnnestMarginalSamples(
    afw::table::Schema const & marginalSchema,
    afw::table::Schema const & directSchema,
    PTR(MarginalSamplingInterpreter) marginalInterpreter,
    PTR(DirectSamplingInterpreter) directInterpreter,
    PTR(afw::math::Random) rng,
    double multiplier
) : _mapper(marginalSchema),
    _multiplier(multiplier),
    _rng(rng),
    _marginalInterpreter(marginalInterpreter),
    _directInterpreter(directInterpreter),
    _prior(_marginalInterpreter->getPrior())
{
    if (_multiplier <= 0.0) {
        _multiplier = static_cast<double>(_directInterpreter->getParameterDim())
            / _marginalInterpreter->getParameterDim();
    }
    if (!_prior) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "Cannot unnest without a Prior"
        );
    }
    MapSamplingFields mapFunc = { &_mapper };
    directSchema.forEach(mapFunc);
    assert(directSchema == _mapper.getOutputSchema()); // should have already thrown an exception if !=
}

void UnnestMarginalSamples::apply(
    ModelFitRecord const & marginalRecord, ModelFitRecord & directRecord
) const {

    // First pass, just to see how many samples we'll have in the new catalog; we want
    // to preallocate space not just for efficiency but so we can call Mixture::updateEM
    // on the columns.
    int nNewSamplesTotal = 0;
    for (
        afw::table::BaseCatalog::const_iterator marginalIter = marginalRecord.getSamples().begin();
        marginalIter != marginalRecord.getSamples().end();
        ++marginalIter
    ) {
        nNewSamplesTotal += boost::math::iround(
            _multiplier
            * marginalRecord.getSamples().size()
            * marginalIter->get(_marginalInterpreter->getWeightKey())
        );
    }

    // Second pass, draw new samples
    int const nonlinearDim = _marginalInterpreter->getNonlinearDim();
    int const amplitudeDim = _marginalInterpreter->getAmplitudeDim();
    int const parameterDim = nonlinearDim + amplitudeDim;
    Vector gradient(amplitudeDim);
    Matrix hessian(amplitudeDim, amplitudeDim);
    directRecord.getSamples().getTable()->preallocate(nNewSamplesTotal);
    ArrayKey amplitudeKey = _directInterpreter->getParameterKey().slice(nonlinearDim, parameterDim);
    for (
        afw::table::BaseCatalog::const_iterator marginalIter = marginalRecord.getSamples().begin();
        marginalIter != marginalRecord.getSamples().end();
        ++marginalIter
    ) {
        Scalar marginalWeight = marginalIter->get(_marginalInterpreter->getWeightKey());
        int nNewSamples = boost::math::iround(
            _multiplier * marginalRecord.getSamples().size() * marginalWeight
        );
        if (nNewSamples == 0) continue;
        Scalar directWeight = marginalWeight / nNewSamples;
        ndarray::Array<Scalar const,1,1> nonlinear
            = marginalIter->get(_marginalInterpreter->getParameterKey());
        _marginalInterpreter->unpackNested(*marginalIter, gradient, hessian);
        ndarray::Array<Scalar,2,2> amplitudes = ndarray::allocate(nNewSamples, amplitudeDim);
        ndarray::Array<Scalar,1,1> weights = ndarray::allocate(nNewSamples);
        weights.deep() = directWeight;
        _prior->drawAmplitudes(gradient, hessian, nonlinear, *_rng, amplitudes, weights, true);
        for (int k = 0; k < nNewSamples; ++k) {
            PTR(afw::table::BaseRecord) newRecord = directRecord.getSamples().addNew();
            newRecord->assign(*marginalIter, _mapper);
            newRecord->set(_directInterpreter->getWeightKey(), weights[k]);
            (*newRecord)[_directInterpreter->getNonlinearKey()] = nonlinear;
            (*newRecord)[amplitudeKey] = amplitudes[k];
        }
    }

    // Now we use the direct interpreter and the new samples to compute a single mean and
    // covariance for the amplitudes; we'll use this to set the mean and covariance of the
    // amplitude parameters for each component in in the Mixture PDF, and then run an E-M
    // update to tweak it up and fill in the covariances between the nonlinear and amplitude
    // parameters.
    PTR(Mixture) marginalPdf = marginalRecord.getPdf();
    ndarray::Array<Scalar,1,1> amplitudeMu
        = _directInterpreter->computeAmplitudeMean(directRecord);
    ndarray::Array<Scalar,2,2> amplitudeSigma
        = _directInterpreter->computeAmplitudeCovariance(directRecord, amplitudeMu);
    Mixture::ComponentList directComponents;
    Vector directMu = Vector::Zero(parameterDim);
    Matrix directSigma = Matrix::Zero(parameterDim, parameterDim);
    for (
        Mixture::iterator marginalIter = marginalPdf->begin();
        marginalIter != marginalPdf->end();
        ++marginalIter
    ) {
        directMu.head(nonlinearDim) = marginalIter->getMu();
        directMu.tail(amplitudeDim) = amplitudeMu.asEigen();
        directSigma.topLeftCorner(nonlinearDim, nonlinearDim) = marginalIter->getSigma();
        directSigma.bottomRightCorner(amplitudeDim, amplitudeDim) = amplitudeSigma.asEigen();
        directComponents.push_back(Mixture::Component(marginalIter->weight, directMu, directSigma));
    }
    PTR(Mixture) directPdf = boost::make_shared<Mixture>(
        parameterDim, boost::ref(directComponents), marginalPdf->getDegreesOfFreedom()
    );
    afw::table::BaseColumnView columns = directRecord.getSamples().getColumnView();
    directPdf->updateEM(columns[_directInterpreter->getParameterKey()],
                        columns[_directInterpreter->getWeightKey()]);
    directRecord.setPdf(directPdf);
}

}}} // namespace lsst::meas::multifit
