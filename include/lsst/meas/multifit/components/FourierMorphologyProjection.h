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
 * Declaration of class FourierMorphologyProjection
 */
#ifndef LSST_MEAS_MULTIFIT_COMPONENTS_FOURIER_MORPHOLOGY_PROJECTION_H
#define LSST_MEAS_MULTIFIT_COMPONENTS_FOURIER_MORPHOLOGY_PROJECTION_H

#include <boost/shared_ptr.hpp>
#include <ndarray/fft.hpp>

#include "lsst/meas/multifit/core.h"
#include "lsst/afw/geom/Extent.h"
#include "lsst/meas/multifit/components/MorphologyProjection.h"
#include "lsst/meas/multifit/components/Morphology.h"

namespace lsst {
namespace meas {
namespace multifit {
namespace components {

/**
 *  A derived MorphologyProjection type designed for use with 
 *  FourierModelProjection.
 */
class FourierMorphologyProjection : public MorphologyProjection {
public:
    typedef boost::shared_ptr<FourierMorphologyProjection> Ptr;
    typedef boost::shared_ptr<FourierMorphologyProjection const> ConstPtr;
    
    /// The padded dimensions of the image representation of the model.
    virtual lsst::afw::geom::Extent2I getDimensions() const = 0;

    /// compute the amount of padding needed for the image representation     
    lsst::afw::geom::Extent2I const getPadding() const {
        lsst::afw::geom::Extent2I kernelDimensions(getKernelDimensions());
        
        return lsst::afw::geom::Extent2I::make(
            kernelDimensions.getX()/2 + 1,
            kernelDimensions.getY()/2 + 1
        ); 
    }

    /**
     *  Compute the derivative of the Fourier-space model with respect to the 
     *  linear parameters.
     */
    virtual ndarray::FourierArray<Pixel,3,3> computeLinearParameterDerivative() = 0;

    /**
     *  Compute the derivative of the Fourier-space model with respect to the 
     *  projected nonlinear parameters.
     */
    virtual ndarray::FourierArray<Pixel,3,3> computeProjectedParameterDerivative() = 0;

protected:

    /**
     * Construct a FourierMorphologyProjection.
     */
    FourierMorphologyProjection(
        Morphology::ConstPtr const & morphology,
        lsst::afw::geom::Extent2I const & kernelDimensions, 
        lsst::afw::geom::AffineTransform const & transform
    ) : MorphologyProjection(morphology,kernelDimensions,transform) {}

};


}}}} // namespace lsst::meas::multifit::components

#endif // !LSST_MEAS_MULTIFIT_COMPONENT_FOURIER_MORPHOLOGY_PROJECTION_H

