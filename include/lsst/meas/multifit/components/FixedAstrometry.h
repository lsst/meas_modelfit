// -*- lsst-c++ -*-
#ifndef LSST_MEAS_MULTIFIT_COMPONENTS_FIXEDASTROMETRY_H
#define LSST_MEAS_MULTIFIT_COMPONENTS_FIXEDASTROMETRY_H

#include <Eigen/Core>
#include <boost/shared_ptr.hpp>

#include "lsst/afw/image/Utils.h"
#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/components/Astrometry.h"

namespace lsst {
namespace meas {
namespace multifit {
namespace components {

class FixedAstrometry : public Astrometry {
public:
    typedef boost::shared_ptr<FixedAstrometry> Ptr;
    typedef boost::shared_ptr<FixedAstrometry const> ConstPtr;
    
    explicit FixedAstrometry(lsst::afw::geom::Point2D const & position) 
      : Astrometry(position) {}

    explicit FixedAstrometry(lsst::afw::coord::Coord const & coord) 
      : Astrometry(coord) {}
    virtual ~FixedAstrometry() {}
    
    virtual DerivativeMatrix const & differentiate() const {
        static DerivativeMatrix i;
        return i;
    }
protected:
    virtual int const getParameterSize() const { return 0; }

    virtual Astrometry::Ptr create(
        boost::shared_ptr<ParameterVector const> const &,
        size_t const &
    ) const {
        return boost::make_shared<FixedAstrometry>(computePosition());
    }
};

}}}} // namespace lsst::meas::multifit::components

#endif // !LSST_MEAS_MULTIFIT_COMPONENTS_FIXEDASTROMETRY_H
