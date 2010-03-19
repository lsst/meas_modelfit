#include "lsst/meas/multifit/SersicCache.h"

#include <cmath>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>

namespace multifit = lsst::meas::multifit;
       
multifit::SersicCache::Ptr multifit::SersicCache::_singleton;

multifit::SersicCache::ConstPtr multifit::SersicCache::getInstance() {
    if(!_singleton) {
        lsst::pex::policy::Policy::Ptr policy(
            lsst::pex::policy::Policy::createPolicy(*getDefaultPolicySource())
        );
        lsst::afw::geom::Extent2D resolution(
            lsst::afw::geom::makeExtentD(
                policy->getDouble("kResolution"), 
                policy->getDouble("sersicIndexResolution")
            )
        );
        lsst::afw::geom::BoxD bounds(
            lsst::afw::geom::makePointD(
                policy->getDouble("kMin"), 
                policy->getDouble("sersicIndexMin")
            ),
            lsst::afw::geom::makePointD(
                policy->getDouble("kMax"),
                policy->getDouble("sersicIndexMax")
            )
        );
        FillFunction::Ptr fillFunction(
            new FillFunction(
                policy->getDouble("epsabs"), 
                policy->getDouble("epsrel"),
                policy->getInt("subintervalLimit")
            )
        );
        _singleton.reset(new SersicCache(bounds, resolution, fillFunction));
    }
    return _singleton;
}
void multifit::SersicCache::FillFunction::IntegralParameters::setN(
    double const & n
) {         
    _n = std::abs(n);
    double twoN = _n * 2;
    _kappa = boost::math::gamma_p_inv(twoN, 0.5);
    _norm = M_PI * twoN * std::pow(_kappa, -twoN) * 
        boost::math::tgamma(twoN);
};

double multifit::SersicCache::FillFunction::sersicFunction(
    double radius, void * parameters
) {
    IntegralParameters const & temp = *static_cast<IntegralParameters*>(
        parameters
    );
   
    gsl_sf_result j0;
    int err = gsl_sf_bessel_J0_e(temp.getK()*radius, &j0);
    if(err || j0.val != j0.val)
        j0.val = 0;

    double exponent = -temp.getKappa() * std::pow(radius, 1.0/temp.getN());
    return (radius*j0.val*std::exp(exponent) / temp.getNorm());
}

double multifit::SersicCache::FillFunction::operator() (double x, double y) const {       
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
    gsl_integration_qagiu(
        &func, 0, _epsabs, _epsrel, _limit, ws,
        &result, &abserr
    );

    gsl_integration_workspace_free(ws);
    gsl_set_error_handler(oldErrorHandler);

    return result;
}
