// -*- lsst-c++ -*-
#ifndef LSST_MEAS_MULTIFIT_COMPONENTS_FIXEDASTROMETRY_H
#define LSST_MEAS_MULTIFIT_COMPONENTS_FIXEDASTROMETRY_H

#include <Eigen/Core>
#include <boost/shared_ptr.hpp>

#include "lsst/afw/image/Utils.h"
#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/Astrometry.h"

namespace lsst {
namespace meas {
namespace multifit {

class ComponentModel;

namespace components {

class FixedAstrometry : public Astrometry {
public:
    typedef boost::shared_ptr<FixedAstrometry> Ptr;
    typedef boost::shared_ptr<FixedAstrometry const> ConstPtr;
    
    explicit FixedAstrometry(Astrometry const &astrometry) :
        Astrometry(astrometry.computePosition()) {}

    explicit FixedAstrometry(lsst::afw::geom::Point2D const & point) :         
        Astrometry(point) {}

    explicit FixedAstrometry(
        boost::shared_ptr<ParameterVector const> const & parameters,
        size_t const & start = 0
    ) : Astrometry(parameters, start, true) { }

    virtual Astrometry::Ptr create const(
        boost::shared_ptr<ParameterVector const> const & parameters,
        size_t const & start
    ) const {
        return boost::make_shared<FixedAstrometry>(parameters, start);
    }

    virtual ~FixedAstrometry() {}

    virtual int const getParameterSize() const { return 0; }

    virtual DerivativeMatrix const & differentiate() const {
        static DerivativeMatrix i;
        return i;
    }
      
    void operator=(FixedAstrometry const & other) { 
        assert(false); }
private:

};

}}}} // namespace lsst::meas::multifit::components

#endif // !LSST_MEAS_MULTIFIT_COMPONENTS_FIXEDASTROMETRY_H
