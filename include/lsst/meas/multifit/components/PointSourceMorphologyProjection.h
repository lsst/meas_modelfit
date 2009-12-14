#ifndef LSST_MEAS_MULTIFIT_COMPONENTS_POINT_SOURCE_MORPHOLOGY_PROJECTION_H
#define LSST_MEAS_MULTIFIT_COMPONENTS_POINT_SOURCE_MORPHOLOGY_PROJECTION_H

#include "lsst/meas/multifit/components/FourierMorphologyProjection.h"

namespace lsst {
namespace meas {
namespace multifit {
namespace components {

class PointSourceMorphology;

class PointSourceMorphologyProjection : public FourierMorphologyProjection {
public:
    typedef boost::shared_ptr<PointSourceMorphologyProjection> Ptr;
    typedef boost::shared_ptr<PointSourceMorphologyProjection const> ConstPtr;

    typedef MorphologyProjection::ParameterJacobianMatrix ParameterJacobianMatrix;
    typedef MorphologyProjection::TransformJacobianMatrix TransformJacobianMatrix;
    
    /// \brief Return the Morphology this is a projection of.
    boost::shared_ptr<PointSourceMorphology const> getMorphology() const { 
        return boost::static_pointer_cast<PointSourceMorphology const>(
            MorphologyProjection::getMorphology()
        );
    }

    /**
     *  \brief Return the matrix that maps the output of
     *  ComponentModelProjection::computeProjectedParameterDerivative()
     *  to the morphology block of the nonlinear parameter derivative.
     */
    virtual ParameterJacobianMatrix const & computeProjectedParameterJacobian() const;

    /**
     *  \brief Return the matrix that deprojects the output of
     *  ComponentModelProjection::computeProjectedParameterDerivative() 
     *  to the morphology terms of the WCS parameter derivative.
     */
    virtual TransformJacobianMatrix const & computeTransformParameterJacobian() const;

    /// \brief Return the padded dimensions of the image representation of the model.
    virtual lsst::afw::geom::Extent2I getDimensions() const {
        int kernelSize = getKernelSize() + getPadding() * 2;
        return lsst::afw::geom::Extent2I(kernelSize, kernelSize);
    }

    /**
     *  Compute the derivative of the Fourier-space model with respect to the linear parameters.
     */
    virtual ndarray::FourierArray<Pixel,3,3> computeLinearParameterDerivative();

    /**
     *  Compute the derivative of the Fourier-space model with respect to the 
     *  projected nonlinear parameters.
     */
    virtual ndarray::FourierArray<Pixel,3,3> computeProjectedParameterDerivative();

    /**
     *  \brief Construct a PointSourceMorphologyProjection.
     */
    PointSourceMorphologyProjection(
        boost::shared_ptr<PointSourceMorphology const> const & morphology,
        int kernelSize, 
        lsst::afw::geom::AffineTransform::ConstPtr const & transform
    );

private:

    friend class PointSourceMorphology;

    ndarray::FourierArray<Pixel,3,3> _linearParameterDerivative;
};

}}}} // namespace lsst::meas::multifit::components

#endif // !LSST_MEAS_MULTIFIT_COMPONENTS_POINT_SOURCE_MORPHOLOGY_PROJECTION_H
