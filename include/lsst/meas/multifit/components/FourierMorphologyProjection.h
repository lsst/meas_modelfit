// -*- lsst-c++ -*-
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
        lsst::afw::geom::AffineTransform::ConstPtr const & transform
    ) : MorphologyProjection(morphology,kernelDimensions,transform) {}

};


}}}} // namespace lsst::meas::multifit::components

#endif // !LSST_MEAS_MULTIFIT_COMPONENT_FOURIER_MORPHOLOGY_PROJECTION_H

