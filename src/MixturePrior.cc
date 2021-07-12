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

#include "Eigen/LU"

#include "ndarray/eigen.h"

#include "lsst/pex/exceptions.h"
#include "lsst/meas/modelfit/MixturePrior.h"
#include "lsst/meas/modelfit/TruncatedGaussian.h"

namespace lsst { namespace meas { namespace modelfit {

//------------- MixturePrior --------------------------------------------------------------------------------

MixturePrior::MixturePrior(std::shared_ptr<Mixture const> mixture, std::string const & tag) :
    Prior(tag), _mixture(mixture)
{}

Scalar MixturePrior::marginalize(
    Vector const & gradient, Matrix const & hessian,
    ndarray::Array<Scalar const,1,1> const & parameters
) const {
    return TruncatedGaussian::fromSeriesParameters(0.0, gradient, hessian).getLogIntegral()
        - std::log(_mixture->evaluate(ndarray::asEigenMatrix(parameters)));
}

Scalar MixturePrior::maximize(
    Vector const & gradient, Matrix const & hessian,
    ndarray::Array<Scalar const,1,1> const & nonlinear,
    ndarray::Array<Scalar,1,1> const & amplitudes
) const {
    TruncatedGaussian tg = TruncatedGaussian::fromSeriesParameters(0.0, gradient, hessian);
    ndarray::asEigenMatrix(amplitudes) = tg.maximize();
    return tg.evaluateLog()(ndarray::asEigenMatrix(amplitudes));
}

Scalar MixturePrior::evaluate(
    ndarray::Array<Scalar const,1,1> const & parameters,
    ndarray::Array<Scalar const,1,1> const & amplitudes
) const {
    if ((ndarray::asEigenArray(amplitudes) < 0.0).any()) {
        return 0.0;
    } else {
        return _mixture->evaluate(ndarray::asEigenMatrix(parameters));
    }
}

void MixturePrior::evaluateDerivatives(
    ndarray::Array<Scalar const,1,1> const & nonlinear,
    ndarray::Array<Scalar const,1,1> const & amplitudes,
    ndarray::Array<Scalar,1,1> const & nonlinearGradient,
    ndarray::Array<Scalar,1,1> const & amplitudeGradient,
    ndarray::Array<Scalar,2,1> const & nonlinearHessian,
    ndarray::Array<Scalar,2,1> const & amplitudeHessian,
    ndarray::Array<Scalar,2,1> const & crossHessian
) const {
    _mixture->evaluateDerivatives(nonlinear, nonlinearGradient, nonlinearHessian);
    amplitudeGradient.deep() = 0.0;
    amplitudeHessian.deep() = 0.0;
    crossHessian.deep() = 0.0;
}

void MixturePrior::drawAmplitudes(
    Vector const & gradient, Matrix const & hessian,
    ndarray::Array<Scalar const,1,1> const & nonlinear,
    afw::math::Random & rng,
    ndarray::Array<Scalar,2,1> const & amplitudes,
    ndarray::Array<Scalar,1,1> const & weights,
    bool multiplyWeights
) const {
    TruncatedGaussian::Sampler sampler
        = TruncatedGaussian::fromSeriesParameters(0.0, gradient, hessian).sample();
    sampler(rng, amplitudes, weights, multiplyWeights);
}

namespace {

class EllipseUpdateRestriction : public Mixture::UpdateRestriction {
public:

    virtual void restrictMu(Vector & mu) const {
        mu[0] = 0.0;
        mu[1] = 0.0;
    }

    virtual void restrictSigma(Matrix & sigma) const {
        sigma(0,0) = sigma(1,1) = 0.5*(sigma(0,0) + sigma(1,1));
        sigma(0,1) = sigma(0,1) = 0.0;
        sigma(0,2) = sigma(2,0) = sigma(1,2) = sigma(2,1) = 0.5*(sigma(0,2) + sigma(1,2));
    }

    EllipseUpdateRestriction() : Mixture::UpdateRestriction(3) {}

};

} // anonymous

Mixture::UpdateRestriction const & MixturePrior::getUpdateRestriction() {
    static EllipseUpdateRestriction const instance;
    return instance;
}

}}} // namespace lsst::meas::modelfit
