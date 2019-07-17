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

#include <limits>

#include "ndarray/eigen.h"

#include "boost/math/special_functions/gamma.hpp"

#include "lsst/afw/geom/Angle.h"
#include "lsst/pex/exceptions.h"
#include "lsst/meas/modelfit/SemiEmpiricalPrior.h"
#include "lsst/meas/modelfit/detail/polynomials.h"
#include "lsst/meas/modelfit/TruncatedGaussian.h"

namespace lsst { namespace meas { namespace modelfit {

namespace {

// Exponential distribution in polar coordinates with a softened core:
// f(x) = A x exp(-sqrt(x^2 + t^2)/s) / s
// with A chosen to integrate to 1 over R2
class SoftenedExponential {
public:

    explicit SoftenedExponential(Scalar sigma, Scalar tau) :
        _sigma(sigma),
        _tau(tau),
        _a(std::exp(tau/sigma)/((sigma + tau)*sigma))
    {}

    Scalar p(Scalar x) const {
        Scalar z = std::sqrt(x*x+_tau*_tau);
        return _a*std::exp(-z/_sigma);
    }

    Scalar d1(Scalar x) const {
        Scalar z = std::sqrt(x*x+_tau*_tau);
        return -_a*std::exp(-z/_sigma)*x/(_sigma*z);
    }

    Scalar d2(Scalar x) const {
        Scalar z = std::sqrt(x*x+_tau*_tau);
        return _a*exp(-z/_sigma)*(x*x*z - _sigma*_tau*_tau)/(_sigma*_sigma*z*z*z);
    }

private:
    Scalar _sigma;
    Scalar _tau;
    Scalar _a;
};

// Ellipticity factor in the (separable) prior; this just translates a prior defined
// in polar coordinates to Cartesian coordinates, with the right differential factors.
class EtaDist {
public:

    explicit EtaDist(SemiEmpiricalPriorControl const & ctrl) :
        _m(ctrl.ellipticitySigma, ctrl.ellipticityCore)
    {}

    Scalar p(Scalar eta1, Scalar eta2) const {
        return _m.p(std::sqrt(eta1*eta1 + eta2*eta2))/(2*geom::PI);
    }

    Eigen::Vector2d d1(Scalar eta1, Scalar eta2) const {
        Scalar eta = std::sqrt(eta1*eta1 + eta2*eta2);
        Eigen::Vector2d r = Eigen::Vector2d::Zero();
        if (eta != 0.0) {
            r[0] = eta1/eta;
            r[1] = eta2/eta;
        }
        return _m.d1(eta)*r/(2*geom::PI);;
    }

    Eigen::Matrix2d d2(Scalar eta1, Scalar eta2) const {
        Scalar eta = std::sqrt(eta1*eta1 + eta2*eta2);
        Eigen::Matrix2d r = Eigen::Matrix2d::Zero();
        if (eta != 0.0) {
            Scalar c1 = eta1/eta;
            Scalar c2 = eta2/eta;
            Scalar d1 = _m.d1(eta)/eta;
            Scalar d2 = _m.d2(eta);
            r(0,0) = c2*c2*d1 + c1*c1*d2;
            r(1,1) = c1*c1*d1 + c2*c2*d2;
            r(0,1) = r(1,0) = c1*c2*(d2 - d1);
        }
        r /= (2*geom::PI);
        return r;
    }

private:
    SoftenedExponential _m;
};

// Cubic polynomial segment defined by constrained values and derivatives at its endpoints.
class Cubic {
public:

    explicit Cubic(Eigen::Matrix<Scalar,4,1,Eigen::DontAlign> const & coeffs) : _coeffs(coeffs) {}

    Scalar p(Scalar x) const {
        return detail::Vandermonde<4>::eval(x).dot(_coeffs);
    }

    Scalar d1(Scalar x) const {
        return detail::Vandermonde<4>::differentiate1(x).dot(_coeffs);
    }

    Scalar d2(Scalar x) const {
        return detail::Vandermonde<4>::differentiate2(x).dot(_coeffs);
    }

    Scalar integrate(Scalar x0, Scalar x1) const {
        return detail::Vandermonde<4>::moment(x0, x1).dot(_coeffs);
    }

private:
    Eigen::Matrix<Scalar,4,1,Eigen::DontAlign> _coeffs;
};

// Trivial constant function
class Constant {
public:

    explicit Constant(Scalar v) : _v(v) {}

    Scalar p(Scalar x) const {
        return _v;
    }

    Scalar d1(Scalar x) const {
        return 0.0;
    }

    Scalar d2(Scalar x) const {
        return 0.0;
    }

    Scalar integrate(Scalar x0, Scalar x1) const {
        return x1 - x0;
    }

private:
    Scalar _v;
};

// Student's T distribution normalized to unit peak amplitude
class StudentsT {
public:

    explicit StudentsT(Scalar mu, Scalar sigma, Scalar nu) : _mu(mu), _sigma(sigma), _nu(nu) {}

    Scalar p(Scalar x) const {
        Scalar z = (x - _mu)/_sigma;
        Scalar a = 1.0 + z*z/_nu;
        return std::pow(a, -0.5*(_nu + 1.0));
    }

    Scalar d1(Scalar x) const {
        Scalar z = (x - _mu)/_sigma;
        Scalar a = 1.0 + z*z/_nu;
        return -z*((_nu + 1.0)/_nu)*std::pow(a, -0.5*(_nu + 3.0))/_sigma;
    }

    Scalar d2(Scalar x) const {
        Scalar z = (x - _mu)/_sigma;
        Scalar a = 1.0 + z*z/_nu;
        return std::pow(a, -0.5*(_nu + 5.0))*(_nu + 1.0)*(a*(_nu + 2.0) - (_nu + 3.0)) / (_nu*_sigma*_sigma);
    }

    Scalar integrate(Scalar x0, Scalar x1) const {
        if (x0 == _mu && x1 == std::numeric_limits<Scalar>::infinity()) {
            return 0.5*std::sqrt(geom::PI*_nu)*boost::math::tgamma_delta_ratio(0.5*_nu, 0.5)*_sigma;
        }
        // Should never get here, since the code in this file will only evaluate at the above endpoints,
        // but we should be careful in case we try to use this elsewhere in the future.
        throw LSST_EXCEPT(pex::exceptions::LogicError, "Integral not implemented");
    }

private:
    Scalar _mu;
    Scalar _sigma;
    Scalar _nu;
};


// Radius factor in the (separable) prior.  That factor is a piecewise function: a Student's T tail
// for large radii, flat for moderate radii, and decreasing quickly to zero as a cubic polynomial
// for very small radii.  This class just combines those piecewise components and delegates to the
// appropriate one for the argument.
struct LogRadiusDist {

    explicit LogRadiusDist(SemiEmpiricalPriorControl const & ctrl) :
        _ramp(detail::solveRampPoly(0.0, 1.0, ctrl.logRadiusMinOuter, ctrl.logRadiusMinInner, 0.0, 0.0)),
        _flat(1.0),
        _tail(ctrl.logRadiusMu, ctrl.logRadiusSigma, ctrl.logRadiusNu),
        _break0(ctrl.logRadiusMinOuter),
        _break1(ctrl.logRadiusMinInner),
        _break2(ctrl.logRadiusMu),
        _norm(
            _ramp.integrate(_break0, _break1) +
            _flat.integrate(_break1, _break2) +
            _tail.integrate(_break2, _break3)
        )
    {}

    Scalar p(Scalar x) const {
        if (x > _break1) {
            if (x < _break2) {
                return _flat.p(x)/_norm;
            } else {
                return _tail.p(x)/_norm;
            }
        } else {
            if (x < _break0) {
                return 0.0;
            } else {
                return _ramp.p(x)/_norm;
            }
        }
    }

    Scalar d1(Scalar x) const {
        if (x > _break1) {
            if (x < _break2) {
                return _flat.d1(x)/_norm;
            } else {
                return _tail.d1(x)/_norm;
            }
        } else {
            if (x < _break0) {
                return 0.0;
            } else {
                return _ramp.d1(x)/_norm;
            }
        }
    }

    Scalar d2(Scalar x) const {
        if (x > _break1) {
            if (x < _break2) {
                return _flat.d2(x)/_norm;
            } else {
                return _tail.d2(x)/_norm;
            }
        } else {
            if (x < _break0) {
                return 0.0;
            } else {
                return _ramp.d2(x)/_norm;
            }
        }
    }

private:

    Cubic _ramp;
    Constant _flat;
    StudentsT _tail;
    Scalar _break0;
    Scalar _break1;
    Scalar _break2;
    static Scalar _break3;
    Scalar _norm;
};

Scalar LogRadiusDist::_break3 = std::numeric_limits<Scalar>::infinity();

} // anonymous

void SemiEmpiricalPriorControl::validate() const {
    if (ellipticitySigma <= 0.0) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            (boost::format("ellipticitySigma must be > 0; got %f") % ellipticitySigma).str()
        );
    }
    if (ellipticityCore <= 0.0) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            (boost::format("ellipticityCore must be > 0; got %f") % ellipticityCore).str()
        );
    }
    if (logRadiusMinInner <= logRadiusMinOuter) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            (boost::format("logRadiusMinInner (%f) must be greater than logRadiusMinOuter (%f)")
             % logRadiusMinInner % logRadiusMinOuter).str()
        );
    }
    if (logRadiusMu <= logRadiusMinInner) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            (boost::format("logRadiusMu (%f) must be greater than logRadiusMinInner (%f)")
             % logRadiusMu % logRadiusMinInner).str()
        );
    }
    if (logRadiusSigma <= 0.0) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            (boost::format("logRadiusSigma must be > 0; got %f") % logRadiusSigma).str()
        );
    }
    if (logRadiusNu <= 0.0) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterError,
            (boost::format("logRadiusNu must be > 0; got %f") % logRadiusNu).str()
        );
    }
}


struct SemiEmpiricalPrior::Impl {

    explicit Impl(SemiEmpiricalPriorControl const & ctrl) : eta(ctrl), lnR(ctrl) {}

    EtaDist eta;
    LogRadiusDist lnR;
};


SemiEmpiricalPrior::SemiEmpiricalPrior(Control const & ctrl) {
    ctrl.validate();
    _impl.reset(new Impl(ctrl));
}


Scalar SemiEmpiricalPrior::evaluate(
    ndarray::Array<Scalar const,1,1> const & nonlinear,
    ndarray::Array<Scalar const,1,1> const & amplitudes
) const {
    if ((ndarray::asEigenArray(amplitudes) < 0.0).any()) {
        return 0.0;
    } else {
        return _impl->eta.p(nonlinear[0], nonlinear[1]) * _impl->lnR.p(nonlinear[2]);
    }
}

void SemiEmpiricalPrior::evaluateDerivatives(
    ndarray::Array<Scalar const,1,1> const & nonlinear,
    ndarray::Array<Scalar const,1,1> const & amplitudes,
    ndarray::Array<Scalar,1,1> const & nonlinearGradient,
    ndarray::Array<Scalar,1,1> const & amplitudeGradient,
    ndarray::Array<Scalar,2,1> const & nonlinearHessian,
    ndarray::Array<Scalar,2,1> const & amplitudeHessian,
    ndarray::Array<Scalar,2,1> const & crossHessian
) const {

    nonlinearGradient.deep() = 0.0;
    nonlinearHessian.deep() = 0.0;
    amplitudeGradient.deep() = 0.0;
    amplitudeHessian.deep() = 0.0;
    crossHessian.deep() = 0.0;

    if ((ndarray::asEigenArray(amplitudes) < 0.0).any()) {
        return;
    }

    // Probability densities of separable components
    Scalar eta0 = _impl->eta.p(nonlinear[0], nonlinear[1]);
    Scalar lnR0 = _impl->lnR.p(nonlinear[2]);

    if (eta0 == 0.0 || lnR0 == 0.0) return;

    // First derivatives of separable components
    Eigen::Vector2d eta1 = _impl->eta.d1(nonlinear[0], nonlinear[1]);
    Scalar lnR1 = _impl->lnR.d1(nonlinear[2]);

    // Second derivatives of separable components
    Eigen::Matrix2d eta2 = _impl->eta.d2(nonlinear[0], nonlinear[1]);
    Scalar lnR2 = _impl->lnR.d2(nonlinear[2]);

    // Fill in full derivatives
    nonlinearGradient[0] = eta1[0]*lnR0;
    nonlinearGradient[1] = eta1[1]*lnR0;
    nonlinearGradient[2] = eta0*lnR1;
    nonlinearHessian(0,0) = eta2(0,0)*lnR0;
    nonlinearHessian(0,1) = eta2(0,1)*lnR0;
    nonlinearHessian(0,2) = eta1[0]*lnR1;
    nonlinearHessian(1,0) = eta2(1,0)*lnR0;
    nonlinearHessian(1,1) = eta2(1,1)*lnR0;
    nonlinearHessian(1,2) = eta1[1]*lnR1;
    nonlinearHessian(2,0) = eta1[0]*lnR1;
    nonlinearHessian(2,1) = eta1[0]*lnR1;
    nonlinearHessian(2,2) = eta0*lnR2;
}

Scalar SemiEmpiricalPrior::marginalize(
    Vector const & gradient, Matrix const & hessian,
    ndarray::Array<Scalar const,1,1> const & nonlinear
) const {
    return TruncatedGaussian::fromSeriesParameters(0.0, gradient, hessian).getLogIntegral()
        - std::log(_impl->eta.p(nonlinear[0], nonlinear[1]) * _impl->lnR.p(nonlinear[2]));
}

Scalar SemiEmpiricalPrior::maximize(
    Vector const & gradient, Matrix const & hessian,
    ndarray::Array<Scalar const,1,1> const & nonlinear,
    ndarray::Array<Scalar,1,1> const & amplitudes
) const {
    TruncatedGaussian tg = TruncatedGaussian::fromSeriesParameters(0.0, gradient, hessian);
    ndarray::asEigenMatrix(amplitudes) = tg.maximize();
    return tg.evaluateLog()(ndarray::asEigenMatrix(amplitudes)) -
        std::log(_impl->eta.p(nonlinear[0], nonlinear[1]) * _impl->lnR.p(nonlinear[2]));
}

void SemiEmpiricalPrior::drawAmplitudes(
    Vector const & gradient, Matrix const & fisher,
    ndarray::Array<Scalar const,1,1> const & nonlinear,
    afw::math::Random & rng,
    ndarray::Array<Scalar,2,1> const & amplitudes,
    ndarray::Array<Scalar,1,1> const & weights,
    bool multiplyWeights
) const {
    throw LSST_EXCEPT(
        pex::exceptions::LogicError,
        "drawAmplitudes not yet implemented for SemiEmpiricalPrior"
    );
}

}}} // namespace lsst::meas::modelfit
