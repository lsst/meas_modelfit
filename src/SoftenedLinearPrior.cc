// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2014 LSST Corporation.
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

#define LSST_MAX_DEBUG 10
#include "lsst/pex/logging/Debug.h"
#include "lsst/pex/exceptions.h"
#include "lsst/meas/multifit/SoftenedLinearPrior.h"
#include "lsst/meas/multifit/TruncatedGaussian.h"

namespace lsst { namespace meas { namespace multifit {

//------------- SoftenedLinearPrior -------------------------------------------------------------------------

// Numerics here are neither as robust as they could be (i.e. we use regular polynomials, not orthogonal ones,
// and don't take pains to keep arguments small), nor as efficient as they could be (not-very-clever integer
// powers, etc), but we don't really need a lot of robustness or efficiency here, so hopefully going for
// clarity here is better.
namespace {

// Class that computes rows of the Vandermonde matrix and related matrices;
// the dot product of these row vectors with the polynomial coefficient
// vectors evaluates the polynomial (or computes a derivative).
template <int N>
struct Vandermonde {

    typedef Eigen::Matrix<double,1,N> RowVector;

    // row vector to evaluate a polynomial
    static RowVector eval(double x) {
        RowVector z = RowVector::Zero();
        double y = 1.0;
        for (int i = 0; i < N; ++i, y *= x) {
            z[i] = y;
        }
        return z;
    }

    // row vector to compute the first derivative of a polynomial
    static RowVector differentiate1(double x) {
        RowVector z = RowVector::Zero();
        double y = 1.0;
        for (int i = 1; i < N; ++i, y *= x) {
            z[i] = i*y;
        }
        return z;
    }

    // row vector to compute the first derivative of a polynomial
    static RowVector differentiate2(double x) {
        RowVector z = RowVector::Zero();
        double y = 1.0;
        for (int i = 2; i < N; ++i, y *= x) {
            z[i] = i*(i-1)*y;
        }
        return z;
    }

    // row vector to compute the integral of p(x) x^m dx from x0 to x1
    static RowVector moment(double x0, double x1, int m=0) {
        RowVector z = RowVector::Zero();
        double y0 = x0;
        double y1 = x1;
        for (int j = 0; j < m; ++j, y0 *= x0, y1 *= x1);
        for (int i = 0; i < N; ++i, y0 *= x0, y1 *= x1) {
            z[i] = (y1 / (i+m+1)) - (y0 / (i+m+1));
        }
        return z;
    }

};

Eigen::Vector4d solveRampPoly(double v0, double v1, double x0, double x1, double s0, double s1) {
    // Solve for the coefficients of a cubic polynomial p(x) that goes from
    // p(x0)=0 to p(x1)=v, with p'(x0)=0 and p'(x1)=s.
    Eigen::Vector4d b;
    Eigen::Matrix4d m;
    // p(x0) = v0
    m.row(0) = Vandermonde<4>::eval(x0);
    b[0] = v0;
    // p(x1) = v1
    m.row(1) = Vandermonde<4>::eval(x1);
    b[1] = v1;
    // p'(x0) = s0
    m.row(2) = Vandermonde<4>::differentiate1(x0);
    b[2] = s0;
    // p'(x1) = s1
    m.row(3) = Vandermonde<4>::differentiate1(x1);
    b[3] = s1;
    return m.fullPivLu().solve(b);
}

} // anonymous

SoftenedLinearPrior::SoftenedLinearPrior(Control const & ctrl) :
    _ctrl(ctrl),
    // we start by making the probability equal to 1 at logRadiusMaxInner; we'll rescale everything later
    _logRadiusP1(ctrl.logRadiusMinMaxRatio),
    _logRadiusSlope((1.0 - ctrl.logRadiusMinMaxRatio) / (ctrl.logRadiusMaxInner - ctrl.logRadiusMinInner)),
    _logRadiusPoly1(
        solveRampPoly(0.0, _logRadiusP1, ctrl.logRadiusMinOuter, ctrl.logRadiusMinInner, 0.0, _logRadiusSlope)
    ),
    _logRadiusPoly2(
        solveRampPoly(1.0, 0.0, ctrl.logRadiusMaxInner, ctrl.logRadiusMaxOuter, _logRadiusSlope, 0.0)
    ),
    _ellipticityPoly(solveRampPoly(1.0, 0.0, ctrl.ellipticityMaxInner, ctrl.ellipticityMaxOuter, 0.0, 0.0))
{
    if (ctrl.logRadiusMinMaxRatio <= 0) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            "logRadiusMinMaxRatio must be > 0"
        );
    }
    if (ctrl.logRadiusMinInner < ctrl.logRadiusMinOuter) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            "logRadiusMinInner must be greater than logRadiusMinOuter"
        );
    }
    if (ctrl.logRadiusMinInner > ctrl.logRadiusMaxInner) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            "logRadiusMinInner must be less than logRadiusMaxInner"
        );
    }
    if (ctrl.logRadiusMaxInner > ctrl.logRadiusMaxOuter) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            "logRadiusMaxInner must be less than logRadiusMaxOuter"
        );
    }
    if (ctrl.ellipticityMaxInner > ctrl.ellipticityMaxOuter) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            "ellipticityMaxInner must be less than ellipticityMaxOuter"
        );
    }
    if (ctrl.ellipticityMaxInner <= 0.0) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            "ellipticityMaxInner must be > 0"
        );
    }
    if (ctrl.ellipticityMaxOuter <= 0.0) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            "ellipticityMaxOuter must be > 0"
        );
    }
    pex::logging::Debug log("meas.multifit.SoftenedLinearPrior");
    log.debug<7>("Constructing new SoftenedLinearPrior");
    // Radius distribution is a trapezoid with cubic-softened sides.  First, the trapezoid:
    double logRadiusCoreIntegral = (ctrl.logRadiusMaxInner - ctrl.logRadiusMinInner)
        * 0.5*(_logRadiusP1 + 1.0);
    // ...the softening cubic on the minimum side:
    double logRadiusMinRampIntegral
        = Vandermonde<4>::moment(ctrl.logRadiusMinOuter, ctrl.logRadiusMinInner).dot(_logRadiusPoly1);
    // ...and the softening cubic on the maximum side:
    double logRadiusMaxRampIntegral
        = Vandermonde<4>::moment(ctrl.logRadiusMaxInner, ctrl.logRadiusMaxOuter).dot(_logRadiusPoly2);
    double logRadiusIntegral = logRadiusCoreIntegral + logRadiusMinRampIntegral + logRadiusMaxRampIntegral;
    _logRadiusMinRampFraction = logRadiusMinRampIntegral / logRadiusIntegral;
    _logRadiusMaxRampFraction = logRadiusMaxRampIntegral / logRadiusIntegral;
    log.debug<7>("logRadiusCoreIntegral=%g, logRadiusMinRampIntegral=%g, logRadiusMaxRampIntegral=%g",
                 logRadiusCoreIntegral, logRadiusMinRampIntegral, logRadiusMaxRampIntegral);

    // ellipticity distribution is a cylinder softened at the outside with a cubic radial profile.
    // First, the cylinder:
    double ellipticityCoreIntegral = 2.0*ctrl.ellipticityMaxInner*ctrl.ellipticityMaxInner;
    // ...and now the softened annulus.  Note that we use moments(..., ..., 1) to compute
    // the integral of [p(e) e de], not just [p(e) de].
    double ellipticityMaxRampIntegral = 2.0*M_PI*Vandermonde<4>::moment(
        ctrl.ellipticityMaxInner,
        ctrl.ellipticityMaxOuter,
        1
    ).dot(_ellipticityPoly);
    double ellipticityIntegral = ellipticityCoreIntegral + ellipticityMaxRampIntegral;
    _ellipticityMaxRampFraction = ellipticityMaxRampIntegral / ellipticityIntegral;
    log.debug<7>("ellipticityCoreIntegral=%g, ellipticityMaxRampIntegral=%g",
                 ellipticityCoreIntegral, ellipticityMaxRampIntegral);

    // now we can compute the full 3-d integral and rescale all the radii probabilities;
    // we'll leave the ellipticity probabilities to max out at 1 since they'll get multiplied
    double fullIntegral = logRadiusIntegral * ellipticityIntegral;
    _logRadiusP1 /= fullIntegral;
    _logRadiusPoly1 /= fullIntegral;
    _logRadiusPoly2 /= fullIntegral;
    _logRadiusSlope /= fullIntegral;
}

Scalar SoftenedLinearPrior::evaluate(
    ndarray::Array<Scalar const,1,1> const & nonlinear,
    ndarray::Array<Scalar const,1,1> const & amplitudes
) const {
    if ((amplitudes.asEigen<Eigen::ArrayXpr>() < 0.0).any()) {
        return 0.0;
    } else {
        return _evaluate(nonlinear);
    }
}

void SoftenedLinearPrior::evaluateDerivatives(
    ndarray::Array<Scalar const,1,1> const & nonlinear,
    ndarray::Array<Scalar const,1,1> const & amplitudes,
    ndarray::Array<Scalar,1,1> const & nonlinearGradient,
    ndarray::Array<Scalar,1,1> const & amplitudeGradient,
    ndarray::Array<Scalar,2,1> const & nonlinearHessian,
    ndarray::Array<Scalar,2,1> const & amplitudeHessian,
    ndarray::Array<Scalar,2,1> const & crossHessian
) const {
    Scalar e1 = nonlinear[0];
    Scalar e2 = nonlinear[1];
    Scalar logRadius = nonlinear[2];
    Scalar ellipticity = std::sqrt(e1*e1 + e2*e2);

    nonlinearGradient.deep() = 0.0;
    nonlinearHessian.deep() = 0.0;
    amplitudeGradient.deep() = 0.0;
    amplitudeHessian.deep() = 0.0;
    crossHessian.deep() = 0.0;

    if ((amplitudes.asEigen<Eigen::ArrayXpr>() < 0.0).any()) {
        return;
    }

    if (logRadius <= _ctrl.logRadiusMinOuter) return;
    if (logRadius >= _ctrl.logRadiusMaxOuter) return;
    if (ellipticity >= _ctrl.ellipticityMaxOuter) return;

    // Starting assumption is that we're in the inner segment in all dimensions
    // (use 'e' and 'r' instead of 'ellipticity', 'logRadius' for brevity, from here on)
    Scalar pr = _logRadiusP1 + logRadius * _logRadiusSlope;
    Scalar pe = 1.0;
    Scalar dpr = _logRadiusSlope;
    Scalar dpe = 0.0;
    Scalar d2pe = 0.0;
    Scalar d2pr = 0.0;

    if (ellipticity > _ctrl.ellipticityMaxInner) { // on ellipticity ramp
        pe = _ellipticityPoly.dot(Vandermonde<4>::eval(ellipticity));
        dpe = _ellipticityPoly.dot(Vandermonde<4>::differentiate1(ellipticity));
        d2pe = _ellipticityPoly.dot(Vandermonde<4>::differentiate2(ellipticity));
    }

    if (logRadius < _ctrl.logRadiusMinInner) { // on logRadius min ramp
        pr = _logRadiusPoly1.dot(Vandermonde<4>::eval(logRadius));
        dpr = _logRadiusPoly1.dot(Vandermonde<4>::differentiate1(logRadius));
        d2pr = _logRadiusPoly1.dot(Vandermonde<4>::differentiate2(logRadius));
    } else if (logRadius > _ctrl.logRadiusMaxInner) { // on logRadius max ramp
        pr = _logRadiusPoly2.dot(Vandermonde<4>::eval(logRadius));
        dpr = _logRadiusPoly2.dot(Vandermonde<4>::differentiate1(logRadius));
        d2pr = _logRadiusPoly2.dot(Vandermonde<4>::differentiate2(logRadius));
    }

    // Make ellipticity nonzero if isn't already, so we can divide by it safely (note that if it is
    // zero, the numerator will always be zero as well).
    ellipticity = std::max(ellipticity, std::numeric_limits<Scalar>::epsilon());

    Scalar de_de1 = e1 / ellipticity;
    Scalar de_de2 = e2 / ellipticity;
    Scalar d2e_de1e1 = de_de2 * de_de2 / ellipticity;
    Scalar d2e_de2e2 = de_de1 * de_de1 / ellipticity;
    Scalar d2e_de1e2 = - de_de1 * de_de2 / ellipticity;

    nonlinearGradient[0] = pr * dpe * de_de1;
    nonlinearGradient[1] = pr * dpe * de_de2;
    nonlinearGradient[2] = pe * dpr;

    nonlinearHessian[0][0] = pr * (d2pe*de_de1*de_de1 + dpe*d2e_de1e1);
    nonlinearHessian[1][1] = pr * (d2pe*de_de2*de_de2 + dpe*d2e_de2e2);
    nonlinearHessian[2][2] = pe * d2pr;
    nonlinearHessian[0][1] = nonlinearHessian[1][0] = pr * (d2pe*de_de1*de_de2 + dpe*d2e_de1e2);
    nonlinearHessian[0][2] = nonlinearHessian[2][0] = dpr * dpe * de_de1;
    nonlinearHessian[1][2] = nonlinearHessian[2][1] = dpr * dpe * de_de2;
}

Scalar SoftenedLinearPrior::marginalize(
    Vector const & gradient, Matrix const & hessian,
    ndarray::Array<Scalar const,1,1> const & nonlinear
) const {
    return TruncatedGaussian::fromSeriesParameters(0.0, gradient, hessian).getLogIntegral()
        - std::log(_evaluate(nonlinear));
}

Scalar SoftenedLinearPrior::maximize(
    Vector const & gradient, Matrix const & hessian,
    ndarray::Array<Scalar const,1,1> const & nonlinear,
    ndarray::Array<Scalar,1,1> const & amplitudes
) const {
    TruncatedGaussian tg = TruncatedGaussian::fromSeriesParameters(0.0, gradient, hessian);
    amplitudes.asEigen() = tg.maximize();
    return tg.evaluateLog()(amplitudes.asEigen()) - std::log(_evaluate(nonlinear));
}

void SoftenedLinearPrior::drawAmplitudes(
    Vector const & gradient, Matrix const & fisher,
    ndarray::Array<Scalar const,1,1> const & nonlinear,
    afw::math::Random & rng,
    ndarray::Array<Scalar,2,1> const & amplitudes,
    ndarray::Array<Scalar,1,1> const & weights,
    bool multiplyWeights
) const {
    // TODO
    throw LSST_EXCEPT(
        pex::exceptions::LogicError,
        "drawAmplitudes not yet implemented for SoftenedLinearPrior"
    );
}

Scalar SoftenedLinearPrior::_evaluate(
    ndarray::Array<Scalar const,1,1> const & nonlinear
) const {

    Scalar e1 = nonlinear[0];
    Scalar e2 = nonlinear[1];
    Scalar logRadius = nonlinear[2];
    Scalar ellipticity = std::sqrt(e1*e1 + e2*e2);

    if (logRadius <= _ctrl.logRadiusMinOuter) return 0.0;
    if (logRadius >= _ctrl.logRadiusMaxOuter) return 0.0;
    if (ellipticity >= _ctrl.ellipticityMaxOuter) return 0.0;

    Scalar p = 0.0;
    if (logRadius < _ctrl.logRadiusMinInner) {
        p = _logRadiusPoly1.dot(Vandermonde<4>::eval(logRadius));
    } else if (logRadius > _ctrl.logRadiusMaxInner) {
        p = _logRadiusPoly2.dot(Vandermonde<4>::eval(logRadius));
    } else {
        p = _logRadiusP1 + logRadius * _logRadiusSlope;
    }

    if (ellipticity > _ctrl.ellipticityMaxInner) {
        p *= _ellipticityPoly.dot(Vandermonde<4>::eval(ellipticity));
    }

    return p;
}

}}} // namespace lsst::meas::multifit
