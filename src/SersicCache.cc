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
 
#include "lsst/meas/multifit/SersicCache.h"

#include <cmath>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>

namespace multifit = lsst::meas::multifit;
      
void multifit::SersicCacheFillFunction::IntegralParameters::setN(
    double const & n
) {         
    _n = std::abs(n);
    double twoN = _n * 2;
    switch (_options.radius) {
    case RADIUS_HALF_INTEGRAL:
        _kappa = boost::math::gamma_p_inv(twoN, 0.5);
        _norm = M_PI * twoN * std::pow(_kappa, -twoN) * 
            boost::math::tgamma(twoN);
        break;
    case RADIUS_HALF_MAX:
        _kappa = M_LN2;
        _norm = 1.0;
        break;
    case RADIUS_NATURAL:
        _kappa = 1.0;
        _norm = 1.0;
        break;
    }
};

double multifit::SersicCacheFillFunction::sersicFunction(
    double radius, void * parameters
) {
    IntegralParameters const & temp = *static_cast<IntegralParameters*>(
        parameters
    );
   
    gsl_sf_result j0;
    int err = gsl_sf_bessel_J0_e(temp.getK()*radius, &j0);
    if(err || j0.val != j0.val)
        j0.val = 0;

    double exponent = temp.getKappa() * (-std::pow(radius, 1.0/temp.getN()));
    double cutoff = 0.0;
    if (temp.getOptions().truncate) {
        cutoff = std::exp(temp.getKappa() * (-std::pow(5.0, 1.0/temp.getN())));
    }
    return (radius*j0.val*((std::exp(exponent) - cutoff)) / temp.getNorm());
}

double multifit::SersicCacheFillFunction::operator() (double x, double y) const {       
    gsl_function func;
    func.function = sersicFunction;
    func.params = static_cast<void*>(&_params);

    // compute parameter dimensions, and allocate grid
    if(y != _lastY) {
        _params.setN(y);
    } 
    _params.setK(x);

    double result, abserr;
    gsl_error_handler_t * oldErrorHandler = gsl_set_error_handler_off();
    
    gsl_integration_workspace * ws = gsl_integration_workspace_alloc(_limit);
    if (_params.getOptions().truncate) {
        gsl_integration_qag(
            &func, 0, _params.getOptions().truncationRadius, _epsabs, _epsrel, _limit,
            GSL_INTEG_GAUSS61, ws, &result, &abserr
        );
    } else {  
        gsl_integration_qagiu(&func, 0, _epsabs, _epsrel, _limit, ws, &result, &abserr);
    }
    gsl_integration_workspace_free(ws);
    gsl_set_error_handler(oldErrorHandler);

    return result;
}

multifit::Cache::ConstPtr multifit::makeSersicCache(lsst::pex::policy::Policy policy) {
    lsst::pex::policy::DefaultPolicyFile defSource(
        "meas_multifit",
        "SersicCacheDict.paf",
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
    SersicCacheFillFunction::Options options;
    std::string radiusString = policy.getString("radius");
    if (radiusString == "HALF_INTEGRAL") {
        options.radius = SersicCacheFillFunction::RADIUS_HALF_INTEGRAL;
    } else if (radiusString == "HALF_MAX") {
        options.radius = SersicCacheFillFunction::RADIUS_HALF_MAX;
    } else if (radiusString == "NATURAL") {
        options.radius = SersicCacheFillFunction::RADIUS_NATURAL;
    } else {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Sersic radius must be one of 'HALF_INTEGRAL', 'HALF_MAX', 'NATURAL'."
        );
    }
    options.truncate = policy.getBool("truncate");
    options.truncationRadius = policy.getDouble("truncationRadius");
    SersicCacheFillFunction fillFunction(
        policy.getDouble("epsabs"), 
        policy.getDouble("epsrel"),
        policy.getInt("subintervalLimit"),
        options
    );

    return Cache::make(
        bounds, resolution, &fillFunction, 
        Cache::FunctorFactory::get(""), 
        Cache::FunctorFactory::get(""),
        "Sersic", false
    );
}
