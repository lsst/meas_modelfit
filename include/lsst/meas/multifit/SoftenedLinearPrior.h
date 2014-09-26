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

#ifndef LSST_MEAS_MULTIFIT_SoftenedLinearPrior_h_INCLUDED
#define LSST_MEAS_MULTIFIT_SoftenedLinearPrior_h_INCLUDED

#include "lsst/pex/config.h"
#include "lsst/meas/multifit/Prior.h"

namespace lsst { namespace meas { namespace multifit {

struct SoftenedLinearPriorControl {

    LSST_CONTROL_FIELD(
        ellipticityMaxOuter, double,
        "Maximum ellipticity magnitude (conformal shear units)"
    );

    LSST_CONTROL_FIELD(
        ellipticityMaxInner, double,
        "Ellipticity magnitude (conformal shear units) at which the softened cutoff begins"
    );

    LSST_CONTROL_FIELD(
        logRadiusMinOuter, double,
        "Minimum ln(radius)"
    );

    LSST_CONTROL_FIELD(
        logRadiusMinInner, double,
        "ln(radius) at which the softened cutoff begins towards the minimum"
    );

    LSST_CONTROL_FIELD(
        logRadiusMaxOuter, double,
        "Maximum ln(radius)"
    );

    LSST_CONTROL_FIELD(
        logRadiusMaxInner, double,
        "ln(radius) at which the softened cutoff begins towards the maximum"
    );

    LSST_CONTROL_FIELD(
        logRadiusMinMaxRatio, double,
        "The ratio P(logRadiusMinInner)/P(logRadiusMaxInner)"
    );

    SoftenedLinearPriorControl() :
        ellipticityMaxOuter(2.001), ellipticityMaxInner(2.0),
        logRadiusMinOuter(-10.001), logRadiusMinInner(-10.0),
        logRadiusMaxOuter(3.001), logRadiusMaxInner(3.0),
        logRadiusMinMaxRatio(1.0)
    {}

};

/**
 *  @brief A prior that's linear in radius and flat in ellipticity, with a cubic roll-off at the edges.
 */
class SoftenedLinearPrior : public Prior {
public:

    typedef SoftenedLinearPriorControl Control;

    explicit SoftenedLinearPrior(Control const & ctrl=Control());

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

    Control const & getControl() const { return _ctrl; }

private:

    Scalar _evaluate(ndarray::Array<Scalar const,1,1> const & nonlinear) const;

    Control _ctrl;
    double _logRadiusP1; // probability value at ln(radius) = ctrl.logRadiusMinInner
    double _logRadiusSlope;
    double _logRadiusMinRampFraction;
    double _logRadiusMaxRampFraction;
    double _ellipticityMaxRampFraction;
    Eigen::Matrix<double,4,1,Eigen::DontAlign> _logRadiusPoly1;
    Eigen::Matrix<double,4,1,Eigen::DontAlign> _logRadiusPoly2;
    Eigen::Matrix<double,4,1,Eigen::DontAlign> _ellipticityPoly;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_SoftenedLinearPrior_h_INCLUDED
