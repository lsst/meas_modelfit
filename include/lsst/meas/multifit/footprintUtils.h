#ifndef LSST_MEAS_MULTIFIT_FOOTPRINT_UTILS_H
#define LSST_MEAS_MULTIFIT_FOOTPRINT_UTILS_H

#include "boost/shared_ptr.hpp"
#include "ndarray.hpp"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/meas/multifit/core.h"

namespace lsst {
namespace meas {
namespace multifit {

lsst::afw::detection::Footprint::Ptr makeFootprint(
    lsst::afw::geom::ellipses::Ellipse const & ellipse
);

template <typename MaskPixel>
lsst::afw::detection::Footprint::Ptr clipAndMaskFootprint(
    lsst::afw::detection::Footprint const & footprint,
    typename lsst::afw::image::Mask<MaskPixel>::Ptr const & mask
);

template <typename ImagePixel, typename MaskPixel, typename VariancePixel>
void compressImage(
    lsst::afw::detection::Footprint const & footprint,
    lsst::afw::image::MaskedImage<ImagePixel, MaskPixel, VariancePixel> const & maskedImage,
    ndarray::Array<Pixel, 1, 1> const & imageDest,
    ndarray::Array<Pixel, 1, 1> const & varianceDest
);

template <typename ImagePixel, typename MaskPixel, typename VariancePixel>
void expandImage(
    lsst::afw::detection::Footprint const & footprint,
    lsst::afw::image::MaskedImage<ImagePixel, MaskPixel, VariancePixel> & maskedImage,
    ndarray::Array<Pixel const, 1, 1> const & imageSrc,
    ndarray::Array<Pixel const, 1, 1> const & varianceSrc
);

}}} //end namespace lsst::meas::multifit

#endif
