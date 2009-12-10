#ifndef LSST_MEAS_MULTIFIT_COMPONENT_MODEL_PROJECTION_H
#define LSST_MEAS_MULTIFIT_COMPONENT_MODEL_PROJECTION_H

#include "lsst/meas/multifit/ComponentModel.h"
#include "lsst/meas/multifit/components/MorphologyProjection.h"
#include "lsst/meas/multifit/matrices.h"

#include "lsst/afw/math/AffineTransform.h"

namespace lsst {
namespace meas {
namespace multifit {

/**
 *  \brief A projection of a ComponentModel to a particular set of observing conditions.
 *  
 *  \sa ComponentModel
 *  \sa ComponentModelFactory
 */
class ComponentModelProjection : 
    public ModelProjection, public lsst::afw::math::Convolvable {       
public:
    typedef boost::shared_ptr<ModelProjection> Ptr;
    typedef boost::shared_ptr<ModelProjection const> ConstPtr;

    static const int WCS_PARAMETER_SIZE = 6;

    /// \brief Return the Model instance this is a projection of.
    ComponentModel::ConstPtr getModel() const {
        return boost::static_pointer_cast<ComponentModel const>(ModelProjection::getModel());
    }

    /// \brief Return the Astrometry object this projection is based on.
    components::Astrometry::ConstPtr getAstrometry() const { 
        return getModel()->getAstrometry(); 
    }

    /// \brief Return the MorphologyProjection object this projection is based on.
    components::MorphologyProjection::ConstPtr getMorphologyProjection() const {
        return _morphologyProjection;
    }

    /// \brief Return the AffineTransform that maps global coordinates to image coordinates.
    lst::afw::math::AffineTransform::ConstPtr getTransform() const { return _transform; }

    /// \brief Return the number of parameters that specify the coordinate transformation.
    virtual int const getWcsParameterSize() const { return WCS_PARAMETER_SIZE; }

    /// \brief Return the number of parameters that specify the PSF.
    virtual int const getPsfParameterSize() const = 0;

protected:

    enum ComponentModelProductFlags {
        TRANSLATION_DERIVATIVE = 1 << 0,
        PROJECTED_PARAMETER_DERIVATIVE = 1 << 1
    };

    /// \brief Construct a projection.
    ComponentModelProjection(
        ComponentModel::ConstPtr const & model,
        Kernel::ConstPtr const & kernel,
        Wcs::ConstPtr const & wcs,
        Footprint::ConstPtr const & footprint,
        double photFactor
    );

    /**
     *  @name LocalizedConvolvableImplementation
     */
    //@{
    virtual lsst::afw::math::PointD getPsfPosition() const { 
        return (*_transform)(getAstrometry()->apply()); 
    }
    //@}

    // ------------------------------------------------------------------- //

    /**
     *  @name ProtectedProductComputers
     *
     *  These are fully implemented by ComponentModel by delegating to _computeTranslationDerivative()
     *  and _computeProjectedParameterDerivative(), and should not generally by reimplemented
     *  by subclasses.
     */
    //@{
    virtual void _computeNonlinearParameterDerivative(ndarray::Array<Pixel,2,2> const & matrix);
    virtual void _computeWcsParameterDerivative(ndarray::Array<Pixel,2,2> const & matrix);
    //@}

    /**
     *  Compute the derivative of the model image with respect to image-coordinate translations
     *  (Footprint-compressed).
     */
    virtual void _computeTranslationDerivative(ndarray::Array<Pixel,2,2> const & matrix) = 0;

    /**
     *  Compute the derivative of the model image with respect to its projected morphology
     *  parameters (Footprint-compressed).
     *
     *  The "projected morphology parameters" correspond to image-frame versions of the
     *  nonlinear morphology parameters, and must transform according to the
     *  MorphologyProjection's getParameterJacobian() and getTransformDerivative() outputs.
     */
    virtual void _computeProjectedParameterDerivative(ndarray::Array<Pixel,2,2> const & matrix) = 0;

    /**
     *  \brief Handle a linear parameter change broadcast from the associated Model.
     *
     *  This propogates the change to the stored MorphologyProjection object.  Subclasses
     *  which override should ensure this implementation is still called.
     */
    virtual void _handleLinearParameterChange() {
        ModelProjection::_handleLinearParameterChange();
        _morphologyProjection->_handleLinearParameterChange();
        _validProducts &= (~TRANSLATION_DERIVATIVE);
        _validProducts &= (~PROJECTED_PARAMETER_DERIVATIVE);
    }

    /**
     *  \brief Handle a nonlinear parameter change broadcast from the associated Model.
     *
     *  This propogates the change to the stored MorphologyProjection object.  Subclasses
     *  which override should ensure this implementation is still called.
     */
    virtual void _handleNonlinearParameterChange() {
        ModelProjection::_handleNonlinearParameterChange();
        _morphologyProjection->_handleMorphologyParameterChange();
        _validProducts &= (~TRANSLATION_DERIVATIVE);
        _validProducts &= (~PROJECTED_PARAMETER_DERIVATIVE);
    }

    components::MorphologyProjection::Ptr _getMorphologyProjection() { 
        return _morphologyProjection;
    }

private:
    
    typedef Eigen::Map< Eigen::Matrix<Pixel,Eigen::Dynamic,2> > TranslationMatrixMap;

    void _ensureTranslationDerivative();
    void _ensureProjectedParameterDerivative();

    MatrixMap getAstrometryParameterMatrixView(ndarray::Array<Pixel,2,2> const & array) {
        return MatrixMap(
            array.getData(), 
            array.getSize<1>(), 
            getModel()->getAstrometry()->getAstrometryParameterSize()
        );
    }

    MatrixMap getMorphologyParameterMatrixView(ndarray::Array<Pixel,2,2> const & array) {
        int offset = getModel()->getAstrometry()->getAstrometryParameterSize() * array.getStride<0>();
        return MatrixMap(
            array.getData() + offset,
            array.getSize<1>(), 
            getModel()->getMorphology()->getMorphologyParameterSize()
        );
    }

    TranslationMatrixMap getTranslationMatrixView() {
        return TranslationMatrixMap(
            _translationDerivative.getData(),
            _translationDerivative.getSize<1>(),
            2
        );
    }

    MatrixMap getProjectedParameterMatrixView() {
        return getMatrixView(_projectedParameterDerivative);
    }

    int _validProducts;
    lsst::afw::math::AffineTransform::ConstPtr _transform; ///< Transform from global coordinates to this projection
    components::MorphologyProjection::Ptr _morphologyProjection;
    ndarray::Array<Pixel,2,2> _translationDerivative;
    ndarray::Array<Pixel,2,2> _projectedParameterDerivative;
};

}}} // namespace multifit

#endif // !LSST_MEAS_MULTIFIT_COMPONENT_MODEL_PROJECTION_H
