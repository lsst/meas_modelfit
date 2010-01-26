// -*- lsst-c++ -*-
#ifndef LSST_MEAS_MULTIFIT_COMPONENTS_MORPHOLOGY_PROJECTION_H
#define LSST_MEAS_MULTIFIT_COMPONENTS_MORPHOLOGY_PROJECTION_H

#include <boost/shared_ptr.hpp>

#include "lsst/meas/multifit/core.h"

namespace lsst {
namespace meas {
namespace multifit {

class ComponentModelProjection;

namespace components {

class Morphology;

/**
 * A projection of a Morphology object.
 *
 * A MorphologyProjection should only exist as a data member of a 
 * ComponentModel.
 */
class MorphologyProjection : private boost::noncopyable {
public:
    typedef boost::shared_ptr<MorphologyProjection> Ptr;
    typedef boost::shared_ptr<MorphologyProjection const> ConstPtr;
    
    typedef Eigen::Matrix<Pixel, Eigen::Dynamic, Eigen::Dynamic> ParameterJacobianMatrix;
    typedef Eigen::Matrix<Pixel,Eigen::Dynamic,6> TransformJacobianMatrix;

    virtual ~MorphologyProjection() {}

    /// Return the Morphology this is a projection of.
    boost::shared_ptr<Morphology const> getMorphology() const { 
        return _morphology; 
    }

    /// Return the kernel dimensions this MorphologyProjection expects.
    lsst::afw::geom::Extent2I const & getKernelDimensions() const { 
        return _kernelDimensions; 
    }

    /**
     * AffineTransform that transforms this projection to the global coordinate 
     * space.
     */
    lsst::afw::geom::AffineTransform::ConstPtr getTransform() const { 
        return _transform; 
    }

    /**
     *  Return the matrix that maps the output of 
     *  ComponentModelProjection::computeProjectedParameterDerivative to the 
     *  morphology block of the nonlinear parameter derivative.
     */
    virtual ParameterJacobianMatrix const & computeProjectedParameterJacobian() const = 0;

    /**
     *  Return the matrix that deprojects the output of 
     *  ComponentModelProjection::computeProjectedParameterDerivative to the 
     *  morphology terms of the WCS parameter derivative.
     */
    virtual TransformJacobianMatrix const & computeTransformParameterJacobian() const = 0;

protected:

    /**
     * Handle a change in the linear parameters, as propogated by the owning
     * ComponentModelProjection.
     */
    virtual void _handleLinearParameterChange() {}

    /**
     * Handle a change in the nonlinear morphology parameters, as propogated 
     * by the owning ComponentModelProjection.
     */
    virtual void _handleNonlinearParameterChange() {}

    /**
     * Construct a MorphologyProjection.
     */
    MorphologyProjection(
        boost::shared_ptr<Morphology const> const & morphology,
        lsst::afw::geom::Extent2I const & kernelDimensions, 
        lsst::afw::geom::AffineTransform::ConstPtr const & transform
    ) : _morphology(morphology), 
        _kernelDimensions(kernelDimensions), 
        _transform(transform) {}

private:
    friend class Morphology;
    friend class lsst::meas::multifit::ComponentModelProjection;

    boost::shared_ptr<Morphology const> _morphology;
    lsst::afw::geom::Extent2I _kernelDimensions;
    lsst::afw::geom::AffineTransform::ConstPtr _transform;
};

}}}} // namespace lsst::meas::multifit::components

#endif // !LSST_MEAS_MULTIFIT_COMPONENTS_MORPHOLOGY_PROJECTION_H

