// -*- lsst-c++ -*-
/**
 * @file
 * Declaration of class PointSourceMorphologyProjection
 */
#ifndef LSST_MEAS_MULTIFIT_COMPONENTS_POINT_SOURCE_MORPHOLOGY_PROJECTION_H
#define LSST_MEAS_MULTIFIT_COMPONENTS_POINT_SOURCE_MORPHOLOGY_PROJECTION_H

#include "lsst/meas/multifit/components/FourierMorphologyProjection.h"

namespace lsst {
namespace meas {
namespace multifit {
namespace components {

class PointSourceMorphology;

/**
 * Subclass of FourierMorphologyProjection used to represent projections of 
 * static point-sources
 */
class PointSourceMorphologyProjection : public FourierMorphologyProjection {
public:
    typedef boost::shared_ptr<PointSourceMorphologyProjection> Ptr;
    typedef boost::shared_ptr<PointSourceMorphologyProjection const> ConstPtr;
    typedef MorphologyProjection::ParameterJacobianMatrix ParameterJacobianMatrix;
    typedef MorphologyProjection::TransformJacobianMatrix TransformJacobianMatrix;
    
    PointSourceMorphologyProjection(
        boost::shared_ptr<PointSourceMorphology const> const & morphology,
        lsst::afw::geom::Extent2I const & kernelDimensions, 
        lsst::afw::geom::AffineTransform::ConstPtr const & transform
    ); 

    /// Imutable access to the PointSourceMorphology this is a projection of.
    boost::shared_ptr<PointSourceMorphology const> getMorphology() const { 
        return boost::static_pointer_cast<PointSourceMorphology const>(
            MorphologyProjection::getMorphology()
        );
    }

    // MorphologyProjection --------------------------------------------------
    virtual ParameterJacobianMatrix const & computeProjectedParameterJacobian() const;
    virtual TransformJacobianMatrix const & computeTransformParameterJacobian() const;

    // FourierMorphologyProjection --------------------------------------------
    virtual lsst::afw::geom::Extent2I getDimensions() const {
        return getKernelDimensions() + getPadding() * 2;
    }

    virtual ndarray::FourierArray<Pixel,3,3> computeLinearParameterDerivative();
    virtual ndarray::FourierArray<Pixel,3,3> computeProjectedParameterDerivative();
    
private:
    friend class PointSourceMorphology;

    ndarray::FourierArray<Pixel,3,3> _linearParameterDerivative;
};

}}}} // namespace lsst::meas::multifit::components

#endif // !LSST_MEAS_MULTIFIT_COMPONENTS_POINT_SOURCE_MORPHOLOGY_PROJECTION_H
