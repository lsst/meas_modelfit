#ifndef LSST_MEAS_MULTIFIT_COMPONENT_MODEL_PROJECTION_H
#define LSST_MEAS_MULTIFIT_COMPONENT_MODEL_PROJECTION_H

#include "lsst/meas/multifit/ComponentModel.h"
#include "lsst/meas/multifit/components/MorphologyProjection.h"
#include "lsst/meas/multifit/matrices.h"

#include "lsst/afw/geom/AffineTransform.h"

namespace lsst {
namespace meas {
namespace multifit {

/**
 *  A projection of a ComponentModel to a particular set of observing conditions.
 *  
 *  @sa ModelProjection
 *  @sa ComponentModel
 *  @sa ComponentModelFactory
 */
class ComponentModelProjection : 
    public ModelProjection  {       
public:
    typedef boost::shared_ptr<ModelProjection> Ptr;
    typedef boost::shared_ptr<ModelProjection const> ConstPtr;

    static const int WCS_PARAMETER_SIZE = 6;

    /// Model instance this projection is associated with.
    ComponentModel::ConstPtr getModel() const {
        return boost::static_pointer_cast<ComponentModel const>(ModelProjection::getModel());
    }

    /// Astrometry object this projection is based on.
    components::Astrometry::ConstPtr getAstrometry() const { 
        return getModel()->getAstrometry(); 
    }

    /**
     * Imutable reference to the MorphologyProjection object this projection is based on.
     */
    components::MorphologyProjection::ConstPtr getMorphologyProjection() const {
        return _morphologyProjection;
    }

    /// AffineTransform that maps global coordinates to image coordinates.
    lsst::afw::geom::AffineTransform::ConstPtr const & getTransform() const { return _transform; }

    /// Number of parameters that specify the coordinate transformation.
    virtual int const getWcsParameterSize() const { return WCS_PARAMETER_SIZE; }

    /// Number of parameters that specify the PSF.
    virtual int const getPsfParameterSize() const = 0;

protected:

    enum ComponentModelProductFlags {
        TRANSLATION_DERIVATIVE = 1 << 0,
        PROJECTED_PARAMETER_DERIVATIVE = 1 << 1
    };

    ComponentModelProjection(
        ComponentModel::ConstPtr const & model,
        PsfConstPtr const & psf,
        WcsConstPtr const & wcs,
        FootprintConstPtr const & footprint
    );

    /**
     * Compute the image-space (xy) coordinates where the PSF should be centered
     */
    virtual lsst::afw::geom::Point2D _getPsfPosition() const { 
        return (*_transform)(getAstrometry()->computePosition()); 
    }

    /**
     *  @name Protected Product Computers
     *
     *  These are fully implemented by ComponentModel by delegating to _computeTranslationDerivative()
     *  and _computeProjectedParameterDerivative(), and should not generally by reimplemented
     *  by subclasses.
     */
    //@{
    virtual void _computeNonlinearParameterDerivative(ndarray::Array<Pixel,2,1> const & matrix);
    virtual void _computeWcsParameterDerivative(ndarray::Array<Pixel,2,1> const & matrix);


    /**
     *  Compute the derivative of the model image with respect to image-coordinate translations
     *  (Footprint-compressed).
     */
    virtual void _computeTranslationDerivative(ndarray::Array<Pixel,2,1> const & matrix) = 0;

    /**
     *  Compute the derivative of the model image with respect to its projected morphology
     *  parameters (Footprint-compressed).
     *
     *  The "projected morphology parameters" correspond to image-frame versions of the
     *  nonlinear morphology parameters, and must transform according to the
     *  MorphologyProjection's getParameterJacobian() and getTransformDerivative() outputs.
     */
    virtual void _computeProjectedParameterDerivative(ndarray::Array<Pixel,2,1> const & matrix) = 0;
    //@}

    /**
     *  @name Product Enabled Checkers
     */
    //@{
    virtual bool hasWcsParameterDerivative() const { return false; }
    virtual bool hasTranslationDerivative() const { return true; }
    virtual bool hasProjectedParameterDerivative() const {
        return getMorphologyProjection()->getMorphology()->getNonlinearParameterSize() > 0;
    }
    //@}
    
    virtual void _handleLinearParameterChange() {
        ModelProjection::_handleLinearParameterChange();
        _morphologyProjection->_handleLinearParameterChange();
    }
    virtual void _handleNonlinearParameterChange() {
        ModelProjection::_handleNonlinearParameterChange();
        _morphologyProjection->_handleNonlinearParameterChange();
    }

    /**
     * Mutable reference to the MorphologyProjection object this projection is based on.
     */
    components::MorphologyProjection::Ptr _getMorphologyProjection() { 
        return _morphologyProjection;
    }

private:
    
    typedef Eigen::Map< Eigen::Matrix<Pixel,Eigen::Dynamic,2> > TranslationMatrixMap;


    void _ensureTranslationDerivative();
    void _ensureProjectedParameterDerivative();

#if 0 // TODO: Reimplement this when Eigen::Map has a constructor that takes strides.
    MatrixMapBlock getAstrometryParameterMatrixView(ndarray::Array<Pixel,2,1> const & array) {
        MatrixMap map(
            array.getData(), 
            array.getStride<0>(), 
            getModel()->getAstrometry()->getAstrometryParameterSize()
        );
        return MatrixMapBlock(map, 0, 0, array.getSize<1>(), map.cols());

    }

    MatrixMapBlock getMorphologyParameterMatrixView(ndarray::Array<Pixel,2,1> const & array) {
        int offset = getModel()->getAstrometry()->getAstrometryParameterSize();
        MatrixMap map(
            array.getData(),
            array.getStride<0>(), 
            getModel()->getMorphology()->getMorphologyParameterSize()
        );
        return MatrixMapBlock(map, 0, offset, array.getSize<1>(), map.cols());
    }

    MatrixMapBlock getProjectedParameterMatrixView() {
        return getMatrixView(_projectedParameterDerivative);
    }
#endif

    TranslationMatrixMap getTranslationMatrixView() {
        return TranslationMatrixMap(
            _translationDerivative.getData(),
            _translationDerivative.getSize<1>(),
            2
        );
    }

    int _validProducts;

    ///< Transform from global coordinates to this projection
    lsst::afw::geom::AffineTransform::ConstPtr _transform; 
    ///< MorphologyProjection this ComponentModelProjection is based on
    components::MorphologyProjection::Ptr _morphologyProjection;

    ndarray::Array<Pixel,2,2> _translationDerivative;
    ndarray::Array<Pixel,2,2> _projectedParameterDerivative;
};

}}} // namespace multifit

#endif // !LSST_MEAS_MULTIFIT_COMPONENT_MODEL_PROJECTION_H
