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
 
#include "lsst/meas/multifit/RobustSersicCache.h"

#include <cmath>
#include <numeric>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/utility.hpp>

#define LSST_MAX_DEBUG 10
#include "lsst/pex/logging/Debug.h"

namespace multifit = lsst::meas::multifit;
namespace pexLog = lsst::pex::logging;

multifit::RobustSersicCacheFillFunction::RobustSersicCacheFillFunction(
    double epsabs, double epsrel, double truncationRadius
) : 
    Cache::FillFunction(0),
    _lastY(std::numeric_limits<double>::quiet_NaN()),
    _epsabs(epsabs), _epsrel(epsrel), _params(truncationRadius)
{
    gsl_sf_result zero;
    gsl_sf_bessel_zero_J0_e(2, &zero);
    _besselZeros.push_back(zero.val);
}

void multifit::RobustSersicCacheFillFunction::IntegralParameters::setN(
    double const & n
) {         
    _n = std::abs(n);
    _exponentA = std::exp(-std::pow(_a, 1.0 / _n));
};

double multifit::RobustSersicCacheFillFunction::integrand(
    double x, void * parameters
) {
    IntegralParameters const & params = *static_cast<IntegralParameters*>(parameters);
   
    gsl_sf_result j0;
    int err = gsl_sf_bessel_J0_e(x, &j0);
    if(err || j0.val != j0.val) {
        pexLog::Debug debug("lsst.meas.multifit.RobustSersicCache");
        debug.debug<3>("Error computing Bessel function at %f: %f", x, j0.val);
        j0.val = 0;
    }
    double exponentX = std::exp(-std::pow(x / params.getK(), 1.0 / params.getN()));
    return x * (exponentX - params.getExponentA()) * j0.val;
}

double multifit::RobustSersicCacheFillFunction::operator() (double x, double y) const {
    pexLog::Debug debug("lsst.meas.multifit.RobustSersicCache");
    gsl_function func;
    func.function = integrand;
    func.params = static_cast<void*>(&_params);

    // compute parameter dimensions
    if(y != _lastY) {
        _params.setN(y);
    } 
    _params.setK(x);

    if (std::fabs(x) < std::numeric_limits<double>::epsilon()) return 0.0;

    // ensure we know the zeros of J_0 up to the integral upper limit
    double max = _params.getK() * _params.getA();
    while (_besselZeros.back() < max) {
        gsl_sf_result zero;
        int err = gsl_sf_bessel_zero_J0_e(2*_besselZeros.size(), &zero);
        if (err) {
            debug.debug<3>("Error computing %dth Bessel function zero: %f",
                           2*_besselZeros.size(), zero.val);
        }
        debug.debug<10>("Extending Bessel function zeros (%dth zero is %f).",
                        2*_besselZeros.size(), zero.val);
        _besselZeros.push_back(zero.val);
    }

    // integrate between each pair of endpoints
    std::vector<double>::iterator lastBesselZero
        = std::lower_bound(_besselZeros.begin(), _besselZeros.end(), max);
    std::vector<double> endpoints((lastBesselZero - _besselZeros.begin()) + 2);
    std::copy(_besselZeros.begin(), lastBesselZero, boost::next(endpoints.begin()));
    endpoints.front() = 0.0;
    endpoints.back() = max;
    std::vector<double> results(endpoints.size() - 1, 0.0);

    double abserr;
    std::size_t nEval;
    gsl_error_handler_t * oldErrorHandler = gsl_set_error_handler_off();

    std::vector<double>::const_iterator i1 = endpoints.begin();
    std::vector<double>::const_iterator i2 = boost::next(i1);
    std::vector<double>::iterator r = results.begin();
    while (i2 != endpoints.end()) {
        if (std::fabs(*i2 - *i1) > std::numeric_limits<double>::epsilon()) {
            gsl_integration_qng(&func, *i1, *i2, _epsabs, _epsrel, &(*r), &abserr, &nEval);
            debug.debug<9>("Integrated from %f to %f with %d function evaluations: %e.", 
                           *i1, *i2, nEval, &(*r));
        }
        ++i1;
        ++i2;
        ++r;
    }
    double total = std::accumulate(results.rbegin(), results.rend(), 0.0);

    debug.debug<9>("Total integral: %e.", total);
    gsl_set_error_handler(oldErrorHandler);

    return total / (_params.getK() * _params.getK());
}

multifit::Cache::ConstPtr multifit::makeRobustSersicCache(lsst::pex::policy::Policy policy) {
    lsst::pex::policy::DefaultPolicyFile defSource(
        "meas_multifit",
        "RobustSersicCacheDict.paf",
        "policy"
    );
    lsst::pex::policy::Policy defPol(defSource);
    policy.mergeDefaults(defPol);

    lsst::afw::geom::Extent2D resolution = lsst::afw::geom::makeExtentD(
        policy.getDouble("kResolution"), 
        policy.getDouble("sersicIndexResolution")
    );
    lsst::afw::geom::BoxD bounds(
        lsst::afw::geom::makePointD(
            policy.getDouble("kMin"), 
            policy.getDouble("sersicIndexMin")
        ),
        lsst::afw::geom::makePointD(
            policy.getDouble("kMax"),
            policy.getDouble("sersicIndexMax")
        )
    );
    RobustSersicCacheFillFunction fillFunction(
        policy.getDouble("epsabs"), 
        policy.getDouble("epsrel"),
        policy.getDouble("truncationRadius")
    );

    return Cache::make(bounds, resolution, &fillFunction, "", "", "RobustSersic", false);
}
