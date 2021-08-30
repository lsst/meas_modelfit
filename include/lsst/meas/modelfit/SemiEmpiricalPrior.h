// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2015-2016 LSST/AURA
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

#ifndef LSST_MEAS_MODELFIT_SemiEmpiricalPrior_h_INCLUDED
#define LSST_MEAS_MODELFIT_SemiEmpiricalPrior_h_INCLUDED

#include "lsst/pex/config.h"
#include "lsst/meas/modelfit/Prior.h"

namespace lsst { namespace meas { namespace modelfit {

struct SemiEmpiricalPriorControl {

    LSST_CONTROL_FIELD(
        ellipticitySigma, double,
        "Width of exponential ellipticity distribution (conformal shear units)."
    );

    LSST_CONTROL_FIELD(
        ellipticityCore, double,
        "Softened core width for ellipticity distribution (conformal shear units)."
    );

    LSST_CONTROL_FIELD(
        logRadiusMinOuter, double,
        "Minimum ln(radius)."
    );

    LSST_CONTROL_FIELD(
        logRadiusMinInner, double,
        "ln(radius) at which the softened cutoff begins towards the minimum"
    );

    LSST_CONTROL_FIELD(
        logRadiusMu, double,
        "Mean of the Student's T distribution used for ln(radius) at large radius, and the transition "
        "point between a flat distribution and the Student's T."
    );

    LSST_CONTROL_FIELD(
        logRadiusSigma, double,
        "Width of the Student's T distribution in ln(radius)."
    );

    LSST_CONTROL_FIELD(
        logRadiusNu, double,
        "Number of degrees of freedom for the Student's T distribution on ln(radius)."
    );

    SemiEmpiricalPriorControl() :
        ellipticitySigma(0.3), ellipticityCore(0.001),
        logRadiusMinOuter(-6.001), logRadiusMinInner(-6.0),
        logRadiusMu(-1.0), logRadiusSigma(0.45), logRadiusNu(50.0)
    {}

    /// Raise InvalidParameterException if the configuration options are invalid.
    void validate() const;

};

/**
 *  @brief A piecewise prior motivated by both real distributions and practical considerations.
 */
class SemiEmpiricalPrior : public Prior {
public:

    typedef SemiEmpiricalPriorControl Control;

    explicit SemiEmpiricalPrior(Control const & ctrl=Control());

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

private:

    struct Impl;

    std::shared_ptr<Impl> _impl;
};

}}} // namespace lsst::meas::modelfit

#endif // !LSST_MEAS_MODELFIT_SemiEmpiricalPrior_h_INCLUDED
