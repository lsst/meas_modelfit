#ifndef LSST_MEAS_MULTIFIT_COMPONENTS_FOURIER_MORPHOLOGY_PROJECTION_H
#define LSST_MEAS_MULTIFIT_COMPONENTS_FOURIER_MORPHOLOGY_PROJECTION_H

#include <boost/shared_ptr.hpp>
#include <ndarray/fft.hpp>
#include <agl/Coordinate.hpp>

#include "multifit/components/MorphologyProjection.hpp"

namespace multifit {
namespace components {

/**
 *  \brief A MorphologyProjection object designed for use with FourierModelProjection.
 */
class FourierMorphologyProjection : public MorphologyProjection {
public:
    typedef boost::shared_ptr<FourierMorphologyProjection> Ptr;
    typedef boost::shared_ptr<FourierMorphologyProjection const> ConstPtr;
    
    /// \brief Return the padded dimensions of the image representation of the model.
    virtual lsst::afw::geom::Extent2I getDimensions() const = 0;

    /**
     *  \brief Return the padding for the image representation (uniform on all sides).
     *
     *  \todo Make sure this calculation is correct.
     */
    int const getPadding() const { return MorphologyProjection::getKernelSize()/2+1; }

    /**
     *  Compute the derivative of the Fourier-space model with respect to the linear parameters.
     */
    virtual ndarray::FourierArray<Pixel,3,3> computeLinearParameterDerivative() = 0;

    /**
     *  Compute the derivative of the Fourier-space model with respect to the 
     *  projected nonlinear parameters.
     */
    virtual ndarray::FourierArray<Pixel,3,3> computeProjectedParameterDerivative() = 0;

protected:

    /**
     *  \brief Construct a MorphologyProjection.
     */
    FourierMorphologyProjection(
        Morphology::ConstPtr const & morphology,
        int kernelSize, 
        lsst::afw::geom::AffineTransform::ConstPtr const & transform
    ) : MorphologyProjection(morphology,kernelSize,transform) {}

};


}}}} // namespace lsst::meas::multifit::components

#endif // !LSST_MEAS_MULTIFIT_COMPONENT_FOURIER_MORPHOLOGY_PROJECTION_H

