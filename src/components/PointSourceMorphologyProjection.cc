// -*- lsst-c++ -*-
/**
 * @file
 * Implementation of class PointSourceMorphologyProjection
 */
#include "lsst/meas/multifit/components/PointSourceMorphology.h"
#include "lsst/meas/multifit/components/PointSourceMorphologyProjection.h"

namespace multifit = lsst::meas::multifit;
namespace components = multifit::components;

components::PointSourceMorphologyProjection::ParameterJacobianMatrixPtr 
components::PointSourceMorphologyProjection::computeProjectedParameterJacobian() const {
    static const ParameterJacobianMatrixPtr m(new ParameterJacobianMatrix()); 
    // matrix has zero size
    return m;
}

components::PointSourceMorphologyProjection::TransformJacobianMatrixPtr
components::PointSourceMorphologyProjection::computeTransformParameterJacobian() const {
    static const TransformJacobianMatrixPtr m(new TransformJacobianMatrix());
    // matrix has zero size
    return m;
}

ndarray::FourierArray<multifit::Pixel,3,3>
components::PointSourceMorphologyProjection::computeLinearParameterDerivative() {
    return _linearParameterDerivative;
}

ndarray::FourierArray<multifit::Pixel,3,3>
components::PointSourceMorphologyProjection::computeProjectedParameterDerivative() {
    return ndarray::FourierArray<Pixel,3,3>();
#if 0
        0,
        ndarray::FourierArray<Pixel,3,3>::Base(
            ndarray::external(
                static_cast<std::complex<Pixel>*>(0),
                ndarray::makeVector(0, getDimensions().getY(), getDimensions().getX())
            )
        )
    );
#endif
}

/**
 * Construct a PointSourceMorphologyProjection.
 */
lsst::meas::multifit::components::PointSourceMorphologyProjection::PointSourceMorphologyProjection(
    PointSourceMorphology::ConstPtr const & morphology,
    lsst::afw::geom::Extent2I const & kernelDimensions, 
    lsst::afw::geom::AffineTransform::ConstPtr const & transform
) : FourierMorphologyProjection(morphology,kernelDimensions,transform),
    _linearParameterDerivative()
{
    ndarray::shallow(_linearParameterDerivative) = ndarray::FourierTransform<Pixel,2>::initializeK(
        ndarray::makeVector(1, getDimensions().getY(), getDimensions().getX())
    );
    _linearParameterDerivative = 1.0;
}
