#ifndef LSST_MEAS_MULTIFIT_TYPEDEFS_H
#define LSST_MEAS_MULTIFIT_TYPEDEFS_H

#include <Eigen/Core>
#include <ndarray_fwd.hpp>
#include <lsst/afw/math/ellipses.hpp>
#include <lsst/afw/image/Exposure.h>
#include <lsst/afw/math/Kernel.h>
#include <lsst/afw/detection/Footprint.h>

namespace lsst {
namespace meas {
namespace multifit {

typedef double Pixel;
typedef double Parameter;
typedef lsst::afw::math::ellipses::LogShear ParameterEllipseCore;
typedef ParameterEllipseCore::Ellipse ParameterEllipse;

typedef Parameter * ParameterIterator;
typedef Parameter const * ParameterConstIterator;
typedef eigen::VectorXd ParameterVector;

typedef ndarray::Array<Pixel, 1, 1> DataVector;
typedef ndarray::Array<Pixel, 2, 2> DataMatrix;

typedef lsst::afw::image::MaskedImage<Pixel> MaskedImage;
typedef lsst::afw::image::Exposure<Pixel> Exposure;
typedef lsst::afw::math::Kernel Kernel;
typedef lsst::afw::detection::Footprint Footprint;

}}} // namespace lsst::meas::multifit

#endif
