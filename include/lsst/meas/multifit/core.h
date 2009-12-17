#ifndef LSST_MEAS_MULTIFIT_CORE_H
#define LSST_MEAS_MULTIFIT_CORE_H

#include <Eigen/Core>

#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/math/Kernel.h"
#include "lsst/afw/detection/Footprint.h"

namespace lsst {
namespace meas {
namespace multifit {

typedef double Pixel;
typedef double Parameter;
typedef lsst::afw::geom::ellipses::LogShear ParameterEllipseCore;
typedef ParameterEllipseCore::Ellipse ParameterEllipse;

typedef Parameter * ParameterIterator;
typedef Parameter const * ParameterConstIterator;
typedef Eigen::VectorXd ParameterVector;
typedef Eigen::aligned_allocator<char> Allocator;

typedef lsst::afw::math::Kernel Kernel;
typedef lsst::afw::image::Wcs Wcs;
typedef lsst::afw::detection::Footprint Footprint;

typedef boost::shared_ptr<Kernel const> KernelConstPtr;
typedef boost::shared_ptr<Wcs const> WcsConstPtr;
typedef boost::shared_ptr<Footprint const> FootprintConstPtr;

}}} // namespace lsst::meas::multifit

#endif
