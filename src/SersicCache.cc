#include "lsst/meas/multifit/SersicCache.h"

#include <cmath>
#include <gsl/gsl_integration.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <Eigen/Core>

namespace multifit = lsst::meas::multifit;

struct IntegralParameters {
public:    
    double const & getK() const {return _k;}
    double const & getN() const {return _n;}
    double const & getKappa() const {return _kappa;}
    double const & getNorm() const {return _norm;}

    void setK(double const & k) {_k = std::abs(k);}
    void setN(double const & n) {        
        _n = std::abs(n);
        double twoN = _n * 2;
        _kappa = boost::math::gamma_p_inv(twoN, 0.5);
        _norm = M_PI * twoN * std::pow(_kappa, -twoN) * boost::math::tgamma(twoN);       
    }
private:
    double _k, _n, _kappa, _norm;
};

double sersicFunction(
    double radius, void * parameters
) {
    IntegralParameters const & temp = *static_cast<IntegralParameters*>(parameters);
    
    double j0 = boost::math::cyl_bessel_j (0, temp.getK()*radius);
    double exponent = -temp.getKappa() * std::pow(radius, 1.0/temp.getN());
    return (radius*j0*std::exp(exponent) / temp.getNorm());
}

void multifit::SersicCache::computeDataPoints() {       
    IntegralParameters params;
    gsl_function func;
    func.function = sersicFunction;
    func.params = static_cast<void*>(&params);

    // compute parameter dimensions, and allocate grid
    int nIndexMax = _dataPoints.rows();
    int kIndexMax = _dataPoints.cols();
    double const & nStep(getRowStep());
    double const & kStep(getColStep());
    
    //declare loop variables
    double n=_parameterBounds.getMinY(), k, result, abserr;
    size_t neval;

    // outer loop over sersic index parameter space (
    for(int i = 0; i < nIndexMax; ++i, n+= nStep) {
        params.setN(n);
        (*_y)[i] = n;

        k = _parameterBounds.getMinX();
        //inner loop over radius parameter space
        for (int j = 0; j < kIndexMax; ++j, k += kStep) {
            params.setK(k);
            (*_x)[j] = k;

            //param 0: function to integrate
            //param 1,2: min, max
            //param 3, 4: absolute and relative epsilon thresholds
            //param 5, 6, 7: output resuls, error, and nEvaluations
            gsl_integration_qng(
                &func,      
                0, GSL_POSINF, 
                _epsabs, _epsrel, 
                &result, &abserr, &neval
            );
           
            //set grid point
            _dataPoints(i,j) = result;
        }
    }
}
       

