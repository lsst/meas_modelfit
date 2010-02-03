#ifndef LSST_MEAS_MULTIFIT_SERSIC_CACHE_H
#define LSST_MEAS_MULTIFIT_SERSIC_CACHE_H

#include "lsst/meas/multifit/Cache.h"

namespace lsst {
namespace meas {
namespace multifit {

class SersicCache : public Cache {
public:
    SersicCache(
        lsst::afw::geom::Extent2I const & dimensions,
        lsst::afw::geom::BoxD const & parameterBounds,
        double const & epsabs, 
        double const & epsrel
    ) : Cache(dimensions, parameterBounds), _epsabs(epsabs), _epsrel(epsrel) {
        computeDataPoints();
    }

protected:
    void computeDataPoints();
private:
    double _epsabs, _epsrel;
};

}}} //namespace lsst::meas::multifit
#endif
