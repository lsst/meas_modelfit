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

#ifndef LSST_MEAS_MULTIFIT_MixturePrior_h_INCLUDED
#define LSST_MEAS_MULTIFIT_MixturePrior_h_INCLUDED

#include "lsst/meas/multifit/Prior.h"
#include "lsst/meas/multifit/Mixture.h"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief A prior that's flat in amplitude parameters, and uses a Mixture for nonlinear parameters.
 */
class MixturePrior : public Prior {
public:

    explicit MixturePrior(PTR(Mixture const) mixture, std::string const & tag="");

    /// @copydoc Prior::evaluate
    virtual Scalar evaluate(
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        ndarray::Array<Scalar const,1,1> const & amplitudes
    ) const;

    /// @copydoc Prior::evaluateDerivatives
    virtual void evaluateDerivatives(
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        ndarray::Array<Scalar const,1,1> const & amplitudes,
        ndarray::Array<Scalar,1,1> const & nonlinearGradient,
        ndarray::Array<Scalar,1,1> const & amplitudeGradient,
        ndarray::Array<Scalar,2,1> const & nonlinearHessian,
        ndarray::Array<Scalar,2,1> const & amplitudeHessian,
        ndarray::Array<Scalar,2,1> const & crossHessian
    ) const;

    /// @copydoc Prior::marginalize
    virtual Scalar marginalize(
        Vector const & gradient, Matrix const & hessian,
        ndarray::Array<Scalar const,1,1> const & nonlinear
    ) const;

    /// @copydoc Prior::maximize
    virtual Scalar maximize(
        Vector const & gradient, Matrix const & hessian,
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        ndarray::Array<Scalar,1,1> const & amplitudes
    ) const;

    /// @copydoc Prior::drawAmplitudes
    virtual void drawAmplitudes(
        Vector const & gradient, Matrix const & fisher,
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        afw::math::Random & rng,
        ndarray::Array<Scalar,2,1> const & amplitudes,
        ndarray::Array<Scalar,1,1> const & weights,
        bool multiplyWeights=false
    ) const;

    /**
     *  @brief Return a MixtureUpdateRestriction appropriate for (e1,e2,r) data.
     *
     *  This restriction object can be used with Mixture<3>::updateEM() to create
     *  a mixture with a strictly isotropic ellipticity distribution.
     */
    static MixtureUpdateRestriction const & getUpdateRestriction();

    PTR(Mixture const) getMixture() const { return _mixture; }

private:
    PTR(Mixture const) _mixture;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_MixturePrior_h_INCLUDED
