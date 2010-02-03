#ifndef LSST_MEAS_MULTIFIT_SERSIC_CACHE_FILL_FUNCTION_H
#define LSST_MEAS_MULTIFIT_SERSIC_CACHE_FILL_FUNCTION_H

#include <limits>
#include <boost/make_shared.hpp>
#include "lsst/meas/multifit/Cache.h"

namespace lsst {
namespace meas {
namespace multifit {


class SersicCacheFillFunction : public Cache::FillFunction {
public:
    SersicCacheFillFunction (
        double const & epsabs, 
        double const & epsrel
    ) : Cache::FillFunction(0),
        _epsabs(epsabs), _epsrel(epsrel), 
        _lastY(std::numeric_limits<double>::quiet_NaN()) 
    {}

    virtual double operator()(double x, double y) const;
    virtual Cache::FillFunction::Ptr clone() const {
        return boost::make_shared<SersicCacheFillFunction>(_epsabs, _epsrel);
    }

private:
    class IntegralParameters {
    public:    
        double const & getK() const {return _k;}
        double const & getN() const {return _n;}
        double const & getKappa() const {return _kappa;}
        double const & getNorm() const {return _norm;}

        void setK(double const & k) {_k = std::abs(k);}
        void setN(double const & n);      
    private:
        double _k, _n, _kappa, _norm;
    };

    static double sersicFunction(double radius, void * parameters);
        
    double _epsabs, _epsrel;
    double _lastY;
    mutable IntegralParameters _params;
};

}}} //namespace lsst::meas::multifit
#endif
