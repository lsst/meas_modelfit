#ifndef LSST_MEAS_MULTIFIT_COMPONENTS_MORPHOLOGY_PROJECTION_H
#define LSST_MEAS_MULTIFIT_COMPONENTS_MORPHOLOGY_PROJECTION_H

#include <boost/shared_ptr.hpp>
#include "lsst/afw/math/ellipses.h"

#include "lsst/meas/multifit/core.h"

namespace meas {
namespace lsst {
namespace multifit {

class ComponentModelProjection;

namespace components {

class Morphology;

/**
 *  \brief A projection of a Morphology object.
 *
 *  A MorphologyProjection should only exist as a data member of a ComponentModel.
 */
class MorphologyProjection : private boost::noncopyable {
public:
    typedef boost::shared_ptr<MorphologyProjection> Ptr;
    typedef boost::shared_ptr<MorphologyProjection const> ConstPtr;

    typedef Eigen::Matrix<Pixel,Eigen::Dynamic,Eigen::Dynamic> ParameterJacobianMatrix;
    typedef Eigen::Matrix<Pixel,Eigen::Dynamic,6> TransformJacobianMatrix;

    virtual ~MorphologyProjection() {}

    /// \brief Return the Morphology this is a projection of.
    boost::shared_ptr<Morphology const> getMorphology() const { return _morphology; }

    /// \brief Return the kernel size this MorphologyProjection expects.
    int const getKernelSize() const { return _kernelSize; }

    /// \brief Return the AffineTransform that relates this projection to the global coordinate frame.
    lsst::afw::math::AffineTransform::ConstPtr getTransform() const { return _transform; }

    /**
     *  \brief Return the matrix that maps the output of
     *  ComponentModelProjection::computeProjectedParameterDerivative()
     *  to the morphology block of the nonlinear parameter derivative.
     */
    virtual ParameterJacobianMatrix const & computeProjectedParameterJacobian() const = 0;

    /**
     *  \brief Return the matrix that deprojects the output of
     *  ComponentModelProjection::computeProjectedParameterDerivative() 
     *  to the morphology terms of the WCS parameter derivative.
     */
    virtual TransformJacobianMatrix const & computeTransformParameterJacobian() const = 0;

protected:

    /**
     *  \brief Handle a change in the linear parameters, as propogated by the owning
     *  ComponentModelProjection.
     */
    virtual void _handleLinearParameterChange() {}

    /**
     *  \brief Handle a change in the (nonlinear) morphology parameters, as propogated by the
     *  owning ComponentModelProjection.
     */
    virtual void _handleMorphologyParameterChange() {}

    /**
     *  \brief Construct a MorphologyProjection.
     */
    MorphologyProjection(
        boost::shared_ptr<Morphology const> const & morphology,
        int kernelSize, 
        lsst::afw::math::AffineTransform::ConstPtr const & transform
    ) : _morphology(morphology), _kernelSize(kernelSize), _transform(transform) {}

private:
    friend class Morphology;
    friend class lsst::meas::multifit::ComponentModelProjection;

    boost::shared_ptr<Morphology const> _morphology;
    int _kernelSize;
    lsst::afw::math::AffineTransform::ConstPtr _transform;
};

}}}} // namespace lsst::meas::multifit::components

#endif // !LSST_MEAS_MULTIFIT_COMPONENTS_MORPHOLOGY_PROJECTION_H

