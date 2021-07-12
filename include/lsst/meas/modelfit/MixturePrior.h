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

#ifndef LSST_MEAS_MODELFIT_MixturePrior_h_INCLUDED
#define LSST_MEAS_MODELFIT_MixturePrior_h_INCLUDED

#include "lsst/meas/modelfit/Prior.h"
#include "lsst/meas/modelfit/Mixture.h"

namespace lsst { namespace meas { namespace modelfit {

/**
 *  @brief A prior that's flat in amplitude parameters, and uses a Mixture for nonlinear parameters.
 */
class MixturePrior : public Prior {
public:

    explicit MixturePrior(std::shared_ptr<Mixture const> mixture, std::string const & tag="");

    /// @copydoc Prior::evaluate
    Scalar evaluate(
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        ndarray::Array<Scalar const,1,1> const & amplitudes
    ) const override;

    /// @copydoc Prior::evaluateDerivatives
    void evaluateDerivatives(
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        ndarray::Array<Scalar const,1,1> const & amplitudes,
        ndarray::Array<Scalar,1,1> const & nonlinearGradient,
        ndarray::Array<Scalar,1,1> const & amplitudeGradient,
        ndarray::Array<Scalar,2,1> const & nonlinearHessian,
        ndarray::Array<Scalar,2,1> const & amplitudeHessian,
        ndarray::Array<Scalar,2,1> const & crossHessian
    ) const override;

    /// @copydoc Prior::marginalize
    Scalar marginalize(
        Vector const & gradient, Matrix const & hessian,
        ndarray::Array<Scalar const,1,1> const & nonlinear
    ) const override;

    /// @copydoc Prior::maximize
    Scalar maximize(
        Vector const & gradient, Matrix const & hessian,
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        ndarray::Array<Scalar,1,1> const & amplitudes
    ) const override;

    /// @copydoc Prior::drawAmplitudes
    void drawAmplitudes(
        Vector const & gradient, Matrix const & fisher,
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        afw::math::Random & rng,
        ndarray::Array<Scalar,2,1> const & amplitudes,
        ndarray::Array<Scalar,1,1> const & weights,
        bool multiplyWeights=false
    ) const override;

    /**
     *  @brief Return a MixtureUpdateRestriction appropriate for (e1,e2,r) data.
     *
     *  This restriction object can be used with Mixture<3>::updateEM() to create
     *  a mixture with a strictly isotropic ellipticity distribution.
     */
    static MixtureUpdateRestriction const & getUpdateRestriction();

    std::shared_ptr<Mixture const> getMixture() const { return _mixture; }

private:
    std::shared_ptr<Mixture const> _mixture;
};

}}} // namespace lsst::meas::modelfit

#endif // !LSST_MEAS_MODELFIT_MixturePrior_h_INCLUDED
