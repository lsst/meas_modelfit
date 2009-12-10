#ifndef LSST_MEAS_MULTIFIT_CORE_H
#define LSST_MEAS_MULTIFIT_CORE_H

#include <Eigen/Core>

#include "lsst/afw/math/ellipses.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/math/Kernel.h"
#include "lsst/afw/detection/Footprint.h"

namespace lsst {
namespace meas {
namespace multifit {

typedef double Pixel;
typedef double Parameter;
typedef lsst::afw::math::ellipses::LogShear ParameterEllipseCore;
typedef ParameterEllipseCore::Ellipse ParameterEllipse;

typedef Parameter * ParameterIterator;
typedef Parameter const * ParameterConstIterator;
typedef Eigen::VectorXd ParameterVector;
typedef Eigen::aligned_allocator<char> Allocator;

typedef lsst::afw::image::MaskedImage<Pixel> MaskedImage;
typedef lsst::afw::image::Exposure<Pixel> Exposure;
typedef lsst::afw::math::Kernel Kernel;
typedef lsst::afw::detection::Footprint Footprint;

}}} // namespace lsst::meas::multifit

#endif
