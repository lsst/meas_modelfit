// -*- lsst-c++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
/**
 * @file
 * Collection of utility functions for creating and working with
 * lsst::afw::detection::Footprint objects
 */
#ifndef LSST_MEAS_MULTIFIT_FOOTPRINT_UTILS_H
#define LSST_MEAS_MULTIFIT_FOOTPRINT_UTILS_H


#include "boost/shared_ptr.hpp"
#include "ndarray.hpp"
#include "lsst/meas/multifit/core.h"

namespace lsst {
namespace meas {
namespace multifit {

lsst::afw::detection::Footprint::Ptr makeFootprint(
    lsst::afw::geom::ellipses::Ellipse const & ellipse
);

template <typename MaskT>
lsst::afw::detection::Footprint::Ptr clipAndMaskFootprint(
    lsst::afw::detection::Footprint const & footprint,
    lsst::afw::image::Mask<MaskT> const & mask,
    MaskT const & bitmask = ~0x0
);

template <typename ImageT, typename MaskT, typename VarianceT>
void compressImage(
    lsst::afw::detection::Footprint const & footprint,
    lsst::afw::image::MaskedImage<ImageT, MaskT, VarianceT> const & maskedImage,
    ndarray::Array<Pixel, 1, 1> const & imageDest,
    ndarray::Array<Pixel, 1, 1> const & varianceDest
);

template <typename ImageT, typename MaskT, typename VarianceT>
void expandImage(
    lsst::afw::detection::Footprint const & footprint,
    lsst::afw::image::MaskedImage<ImageT, MaskT, VarianceT> & maskedImage,
    ndarray::Array<Pixel const, 1, 1> const & imageSrc,
    ndarray::Array<Pixel const, 1, 1> const & varianceSrc
);

template <typename ImageT>
void compressImage(
    lsst::afw::detection::Footprint const & footprint,
    lsst::afw::image::Image<ImageT> const & image,
    ndarray::Array<Pixel, 1, 1> const & imageDest
);

template <typename ImageT>
void expandImage(
    lsst::afw::detection::Footprint const & footprint,
    lsst::afw::image::Image<ImageT> & image,
    ndarray::Array<Pixel const, 1, 1> const & imageSrc
);

}}} //end namespace lsst::meas::multifit

#endif
