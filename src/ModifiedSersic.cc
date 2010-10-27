// -*- lsst-c++ -*-
/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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

#include "lsst/meas/multifit/ModifiedSersic.h"

#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/utility.hpp>
#include <Eigen/Core>
#include <Eigen/Array>
#include <Eigen/LU>

#define LSST_MAX_DEBUG 10
#include "lsst/pex/logging/Debug.h"

namespace multifit = lsst::meas::multifit;
namespace pexLog = lsst::pex::logging;

namespace {

static void gslHandleError(char const * reason, char const * file, int line, int gsl_errno) {
    pexLog::Debug debug("lsst.meas.multifit.ModifiedSersic");
    debug.debug<5>("GSL error encountered: %s in file '%s', line %d.", reason, file, line);
}

/**
 *  A 1d polynomial or rational function model with fixed degree.
 */
template <int N>
class RationalModel {
public:
    typedef Eigen::Matrix<double,N,1> Vector;
    typedef Eigen::Matrix<double,N,N> Matrix;

    /**
     *  @brief Construct a model from the given x,y data points.
     *
     *  @param[in] degree    degree of the highest-order negative exponent term; for example,
     *                       for (a*x + b + c/x + d/x^2), N==4 and degree==2. 
     */
    RationalModel(Vector const & x, Vector const & y, int degree=(N-1)) : _degree(degree) {
        Matrix vandermonde;
        vandermonde.col(N-1).fill(1.0);
        for (int n = N-2; n >= 0; --n)
            vandermonde.col(n) = vandermonde.col(n+1).cwise() * x;
        Vector q(y.cwise() * (x.cwise().pow(_degree)));
        Eigen::LU<Matrix> lu(vandermonde);
        lu.solve(q, &_model);
    }

    double operator()(double x) const {
        Vector p;
        p[N-1] = std::pow(x, -_degree);
        for (int n = N-2; n >= 0; --n)
            p[n] = p[n+1] * x;
        return p.dot(_model);
    }

private:
    int _degree;
    Vector _model;
};

/**
 *  Class that actually does the work for ModifiedSersicHankelTransform::operator().
 *
 *  A new instance is constructed for every different k value.  The constructor
 *  does all the work.
 */
class Integrator {
public:

    static int const INITIAL = 20;
    static int const MAX_TINY_STEPS = 3;

    Integrator(
        multifit::ModifiedSersicFunction const & function,
        std::vector<double> & endpoints,
        std::vector<double> & values,
        double epsrel,
        double epabs,
        double k,
        bool doInterpolation = false
    );
    
    double getFinalTolerance() const { return _epsabs; }

private:

    /// Function to be passed to GSL integrator.
    static double integrand(double x, void * p);

    /// Return the zeros of the zeroth-order Bessel function of the first kind (J_0).
    static std::vector<double> const & getBesselZeros(double max);

    int setupVectors(double max);

    /*
     *  Between iLower and iUpper, use a rational function to estimate the value at (iLower + iUpper)/2.
     *  Check the estimate by actually integrating at that point; if it's close enough, use a higher-order
     *  rational function to fill in the values.  If it's not close enough, subdivide and recurse.
     */
    void recurse(int iLower, double rLower, int iUpper, double rUpper);
    
    /**
     *  Integrate from _endpoints[i] to _endpoints[i+1], assigning to the result to _integrals[i].
     *
     *  Returns the estimated absolute error in the integral, NOT the value
     *  of the integral.
     */
    double integrate(int i);

    multifit::ModifiedSersicFunction const & _function;
    std::vector<double> & _endpoints;
    std::vector<double> & _integrals;
    int _size;
    pexLog::Debug _debug;
    double _epsrel;
    double _epsabs;
    double _k;
    double _kSquared;
};

Integrator::Integrator(
    multifit::ModifiedSersicFunction const & function,
    std::vector<double> & endpoints,
    std::vector<double> & integrals,
    double epsrel,
    double epsabs1,
    double k,
    bool doInterpolation
) : _function(function), _endpoints(endpoints), _integrals(integrals),
    _debug("lsst.meas.multifit.ModifiedSersic"),
    _epsrel(epsrel), _epsabs(epsabs1), _k(k), _kSquared(k*k)
{
    gsl_error_handler_t * oldHandler = gsl_set_error_handler(&gslHandleError);
    double max = _k * _function.getOuter();
    _size = setupVectors(max);
    // Explicitly integrate the first few values, increasing epsabs to match best possible final precision.
    int iLower = 0;
    int tinySteps = 0;
    double initialValue = 0.0;
    double initialVariance = 0.0;
    double epsabs2 = 0.0;
    for (int const stop = std::min(INITIAL, _size); iLower < stop; ++iLower) {
        double sigma = integrate(iLower);
        initialValue += _integrals[iLower];
        initialVariance += sigma;
        epsabs2 = std::max(
            std::max(epsabs2, std::sqrt(initialVariance) / _size),
            std::numeric_limits<double>::epsilon() * std::fabs(initialValue)
        );
        _epsabs = std::max(epsabs1, epsabs2);
        if (std::fabs(initialValue) < _epsabs && std::fabs(_integrals[iLower]) < _epsabs) {
            if (++tinySteps > MAX_TINY_STEPS) {
                _debug.debug<7>("Integral clamped to zero at section %d: current value is %e, "
                                "last term is %e, epsabs is %e.", 
                                iLower, initialValue, _integrals[iLower], _epsabs);
                gsl_set_error_handler(oldHandler);
                return;
            }
        }
    }
    if (!doInterpolation) {
        for (; iLower < _size; ++iLower) {
            integrate(iLower);
        }
    }
    if (iLower == _size) {
        gsl_set_error_handler(oldHandler);
        return;
    }
    // Last value is special; it's between a Bessel zero and max, not a pair of Bessel zeros.
    int iUpper = _size - 1;
    integrate(iUpper);
    if (iLower == iUpper) {
        gsl_set_error_handler(oldHandler);
        return;
    }
    integrate(--iUpper);
    --iLower;
    recurse(
        iLower, 0.5 * (_endpoints[iLower] + _endpoints[iLower + 1]), 
        iUpper, 0.5 * (_endpoints[iUpper] + _endpoints[iUpper + 1])
    );
}

double Integrator::integrand(double x, void * p) {
    Integrator const & self = *reinterpret_cast<Integrator*>(p);
    return self._function(x / self._k) * x * gsl_sf_bessel_J0(x);
}

std::vector<double> const & Integrator::getBesselZeros(double max) {
    static std::vector<double> besselZeros(1, gsl_sf_bessel_zero_J0(1));
    while (besselZeros.back() < max) {
        besselZeros.push_back(gsl_sf_bessel_zero_J0(besselZeros.size() + 1));
    }
    return besselZeros;
}

int Integrator::setupVectors(double max) {
    std::vector<double> const & besselZeros = getBesselZeros(max);
    _endpoints.clear();
    _endpoints.push_back(0.0);
    std::vector<double>::const_iterator i = besselZeros.begin();
    if (*i >= max) {
        _endpoints.push_back(max);
    } else {
        while (true) {
            _endpoints.push_back(*i);
            if (*(++i) >= max) {
                _endpoints.back() = max;
                break;
            }
            if (*(++i) >= max) {
                _endpoints.push_back(max);
                break;
            }
        }
    }
    _integrals.resize(_endpoints.size() - 1);
    std::fill(_integrals.begin(), _integrals.end(), 0.0);
    return _integrals.size();
}

void Integrator::recurse(int iLower, double rLower, int iUpper, double rUpper) {
    if (iUpper - iLower <= 1) return;
    if (iUpper - iLower == 2) {
        integrate(iLower+1);
        return;
    }
    double iMid = (iLower + iUpper) / 2;
    double rMid = 0.5 * (_endpoints[iMid], _endpoints[iMid+1]);
    RationalModel<2> model2(
        RationalModel<2>::Vector(rLower, rUpper),
        RationalModel<2>::Vector(_integrals[iLower], _integrals[iUpper])
    );
    double estimate = model2(rMid);
    double integrationError = integrate(iMid);
    double absError = std::fabs(estimate - _integrals[iMid]);
    double relError = absError / std::fabs(_integrals[iMid]);
    if (absError < _epsabs || relError < _epsrel || absError < integrationError) {
        RationalModel<3> model3(
            RationalModel<3>::Vector(rLower, rMid, rUpper),
            RationalModel<3>::Vector(_integrals[iLower], _integrals[iMid], _integrals[iUpper])
        );
        for (int i = iLower + 1; i != iMid; ++i) {
            _integrals[i] = model3(0.5*(_endpoints[i] + _endpoints[i+1]));
        }
        for (int i = iMid + 1; i != iUpper; ++i) {
            _integrals[i] = model3(0.5*(_endpoints[i] + _endpoints[i+1]));
        }        
    } else {
        recurse(iLower, rLower, iMid, rMid);
        recurse(iMid, rMid, iUpper, rUpper);
    }
}

double Integrator::integrate(int i) {
    gsl_function gslIntegrand;
    gslIntegrand.function = &integrand;
    gslIntegrand.params = this;
    double abserr;
    std::size_t neval;
    int err = gsl_integration_qng(
        &gslIntegrand, _endpoints[i], _endpoints[i+1], 
        _kSquared*_epsabs, _epsrel, 
        &_integrals[i], &abserr, &neval
    );
    abserr /= _kSquared;
    _integrals[i] /= _kSquared;
    if (err) {
        _debug.debug<6>("Error integrating from %f to %f; result=%g, epsabs=%g, epsrel=%g.",
                        _endpoints[i], _endpoints[i+1], _integrals[i], _epsabs, _epsrel);
    }
    return abserr;
}

int const Integrator::INITIAL;
int const Integrator::MAX_TINY_STEPS;

} // unnamed namespace

multifit::ModifiedSersicFunction::ModifiedSersicFunction(double n, double inner, double outer) 
    :  _inner(inner), _outer(outer) 
{
    setSersicIndex(n);
}

void multifit::ModifiedSersicFunction::setSersicIndex(double n) {
    _n = n;
    double e = -std::pow(_inner, 1.0 / _n);
    double f = std::exp(e);
    double df = e * f / (_n * _inner);
    double dr = _inner - _outer;
    _a[1] = 0.5 * df / dr;
    _a[0] = f - _a[1] * dr * dr;
    _a[2] = -2.0 * _a[1] * _outer;
    _a[3] = _a[1] * _outer * _outer;
}


double multifit::ModifiedSersicFunction::operator()(double r) const {
    return (r <= _inner) ? 
        (std::exp(-std::pow(r, 1.0 / _n)) - _a[0])
        : ((r > _outer) ? 0.0 : ((_a[1] * r + _a[2]) * r + _a[3]));
}

double multifit::ModifiedSersicFunction::integrate(double radius, int m) const {
    if (radius < _inner) return integrateInner(radius, m);
    if (radius < _outer) return integrateInner(_inner, m) + integrateOuter(radius, m);
    return integrateInner(_inner, m) + integrateOuter(_outer, m);
}

double multifit::ModifiedSersicFunction::integrateInner(double radius, int m) const {
    gsl_error_handler_t * oldHandler = gsl_set_error_handler(&gslHandleError);
    double ln_gamma = gsl_sf_lngamma(_n * (m+1));
    double ln_inc_gamma = std::log(gsl_sf_gamma_inc_P(_n * (m+1), std::pow(radius, 1.0 / _n)));
    gsl_set_error_handler(oldHandler);
    return _n * std::exp(ln_gamma + ln_inc_gamma) - _a[0] * radius;
}

double multifit::ModifiedSersicFunction::integrateOuter(double radius, int m) const {
    return _a[1] * (std::pow(radius, m+3) / (m+3) - std::pow(_inner, m+3) / (m+3))
        +  _a[2] * (std::pow(radius, m+2) / (m+2) - std::pow(_inner, m+2) / (m+2))
        +  _a[3] * (std::pow(radius, m+1) / (m+1) - std::pow(_inner, m+1) / (m+1));            
}

multifit::ModifiedSersicHankelTransform::ModifiedSersicHankelTransform(
    ModifiedSersicFunction const & function,
    double epsrel, double epsabs, bool doInterpolation
) : _function(function), _epsrel(epsrel), _epsabs(epsabs), _doInterpolation(doInterpolation) {
    if (_epsrel < 0.0) _epsrel = std::sqrt(std::numeric_limits<double>::epsilon());
    if (_epsabs < 0.0) _epsabs = std::numeric_limits<double>::epsilon();
}

double multifit::ModifiedSersicHankelTransform::operator()(double k) const {
    if (k <= _epsabs) {
        return _function.integrate(_function.getOuter(), 1) * gsl_sf_bessel_J0(0.0);
    }
    Integrator integrator(_function, _endpoints, _integrals, _epsrel, _epsabs, k, _doInterpolation);
    double result = sum(_integrals);
    if (std::fabs(result) < integrator.getFinalTolerance())
        result = 0.0;
    return result;
}

double multifit::ModifiedSersicHankelTransform::sum(std::vector<double> & v) {
    std::sort(v.begin(), v.end());
    return std::accumulate(v.begin(), v.end(), static_cast<long double>(0.0));
}
