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

#include "lsst/meas/multifit/DirectSamplingInterpreter.h"
#include "lsst/meas/multifit/ModelFitRecord.h"

namespace lsst { namespace meas { namespace multifit {

namespace {

class DirectSamplingObjective : public SamplingObjective {
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
        if (getInterpreter()->getPrior()) {
            chiSq -= std::log(getInterpreter()->getPrior()->evaluate(nonlinear, amplitudes));
        }
        return chiSq;
    }

    DirectSamplingObjective(
        PTR(SamplingInterpreter) interpreter,
        PTR(Likelihood) likelihood
    ) :
        SamplingObjective(interpreter, likelihood)
    {}

private:
    ndarray::Array<Pixel,1,1> _residuals;
};

Model::NameVector concatenateNameVectors(Model::NameVector const & a, Model::NameVector const & b) {
    Model::NameVector r;
    r.reserve(a.size() + b.size());
    r.insert(r.end(), a.begin(), a.end());
    r.insert(r.end(), b.begin(), b.end());
    return r;
}

} // anonymous

DirectSamplingInterpreter::DirectSamplingInterpreter(
    afw::table::Schema & sampleSchema,
    PTR(Model) model,
    PTR(Prior) prior
) :
    SamplingInterpreter(
        sampleSchema,
        concatenateNameVectors(model->getNonlinearNames(), model->getAmplitudeNames()),
        model,
        prior
    )
{
    _parameterKey = sampleSchema.addField(
        afw::table::Field< afw::table::Array<Scalar> >(
            "parameters", "nonlinear and amplitude parameters at this sample point",
            model->getNonlinearDim() + model->getAmplitudeDim()
        ),
        true // doReplace
    );
    _nonlinearKey = _parameterKey.slice(0, model->getNonlinearDim());
}

ndarray::Array<Scalar,1,1> DirectSamplingInterpreter::computeAmplitudeQuantiles(
    ModelFitRecord const & record,
    ndarray::Array<Scalar const,1,1> const & fractions,
    int index
) const {
    return computeSampleQuantiles(
        record.getSamples(),
        fractions,
        _parameterKey[index + getNonlinearDim()]
    );
}

ndarray::Array<Scalar,1,1> DirectSamplingInterpreter::computeAmplitudeMean(
    ModelFitRecord const & record
) const {
    return computeSampleMean(
        record.getSamples(),
        _parameterKey.slice(getNonlinearDim(), getParameterDim())
    );
}

ndarray::Array<Scalar,2,2> DirectSamplingInterpreter::computeAmplitudeCovariance(
    ModelFitRecord const & record,
    ndarray::Array<Scalar const,1,1> const & mean
) const {
    return computeSampleCovariance(
        record.getSamples(),
        mean,
        _parameterKey.slice(getNonlinearDim(), getParameterDim())
    );
}

PTR(SamplingObjective) DirectSamplingInterpreter::makeObjective(
    PTR(SamplingInterpreter) self,
    PTR(Likelihood) likelihood
) const {
    return boost::make_shared<DirectSamplingObjective>(self, likelihood);
}

void DirectSamplingInterpreter::_packParameters(
    ndarray::Array<Scalar const,1,1> const & nonlinear,
    ndarray::Array<Scalar const,1,1> const & amplitudes,
    ndarray::Array<Scalar,1,1> const & parameters
) const {
    parameters[ndarray::view(0, getNonlinearDim())] = nonlinear;
    parameters[ndarray::view(getNonlinearDim(), getParameterDim())] = amplitudes;
}

void DirectSamplingInterpreter::_unpackNonlinear(
    ndarray::Array<Scalar const,1,1> const & parameters,
    ndarray::Array<Scalar,1,1> const & nonlinear
) const {
    nonlinear.deep() = parameters[ndarray::view(0, getNonlinearDim())];
}

}}} // namespace lsst::meas::multifit
