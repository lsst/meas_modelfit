#include "lsst/meas/multifit/SersicCacheFillFunction.h"

#include <cmath>
#include <gsl/gsl_integration.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>

namespace multifit = lsst::meas::multifit;
        
void multifit::SersicCacheFillFunction::IntegralParameters::setN(
    double const & n
) {         
    _n = std::abs(n);
    double twoN = _n * 2;
    _kappa = boost::math::gamma_p_inv(twoN, 0.5);
    _norm = M_PI * twoN * std::pow(_kappa, -twoN) * 
        boost::math::tgamma(twoN);
};

double multifit::SersicCacheFillFunction::sersicFunction(
    double radius, void * parameters
) {
    IntegralParameters const & temp = *static_cast<IntegralParameters*>(
        parameters
    );
    
    double j0 = boost::math::cyl_bessel_j (0, temp.getK()*radius);
    double exponent = -temp.getKappa() * std::pow(radius, 1.0/temp.getN());
    return (radius*j0*std::exp(exponent) / temp.getNorm());
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
    size_t neval;

    gsl_integration_qng(
        &func,      
        0, GSL_POSINF, 
        _epsabs, _epsrel, 
        &result, &abserr, &neval
    );

    return result;
}
