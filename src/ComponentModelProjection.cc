#include "lsst/meas/multifit/ComponentModelProjection.h"
#include "lsst/meas/multifit/matrices.h"

namespace multifit = lsst::meas::multifit;

/**
 * Construct a ComponentModelProjection
 *
 */
multifit::ComponentModelProjection::ComponentModelProjection(
    ComponentModel::ConstPtr const & model,
    PsfConstPtr const & psf,
    WcsConstPtr const & wcs,
    FootprintConstPtr const & footprint
) : ModelProjection(model, wcs, footprint),
    _validProducts(0),
    _morphologyProjection(), 
    _translationDerivative(), 
    _projectedParameterDerivative()
{
    lsst::afw::geom::PointD center(getAstrometry()->computePosition());
    lsst::afw::geom::AffineTransform wcsTransform(wcs->linearizeAt(center));
    
    _transform = boost::make_shared<lsst::afw::geom::AffineTransform>(
        wcsTransform.invert()
    );

    _morphologyProjection = model->getMorphology()->makeProjection(
        lsst::afw::geom::ExtentI::make(
            psf->getWidth(), psf->getHeight()
        ), 
        _transform
    );
}

void multifit::ComponentModelProjection::_computeNonlinearParameterDerivative(
    ndarray::Array<Pixel,2,1> const & matrix
) {
    if (hasTranslationDerivative()) {
        _ensureTranslationDerivative();
        // TODO: Move this into an inline function when possible.
        MatrixMap astrometryMap(
            matrix.getData(),
            matrix.getStride<0>(),
            getModel()->getAstrometry()->getParameterSize()
        );                      
        MatrixMapBlock astrometryView(
            astrometryMap, 
            0, 0, 
            matrix.getSize<1>(), matrix.getSize<0>()
        );
        // END TODO
        TranslationMatrixMap translationView(getTranslationMatrixView());
        components::Astrometry::DerivativeMatrix astrometryDerivative(
            getModel()->getAstrometry()->differentiate()
        );
        astrometryView = translationView * 
            _transform->getEigenTransform().linear() * astrometryDerivative;
    }
    if (hasProjectedParameterDerivative()) {
        _ensureProjectedParameterDerivative();
        // TODO: Move this into an inline function when possible.
        int offset = getModel()->getAstrometry()->getParameterSize();
        MatrixMap morphologyMap(
            matrix.getData(),
            matrix.getStride<0>(), 
            getModel()->getMorphology()->getNonlinearParameterSize()
        );
        MatrixMapBlock morphologyView(
            morphologyMap, 
            0, offset, 
            matrix.getSize<1>(), morphologyMap.cols()
        );
        // END TODO
        // TODO: Move this into an inline function when possible.
        MatrixMap projectedMap(
            _projectedParameterDerivative.getData(),
            _projectedParameterDerivative.getStride<0>(), 
            _projectedParameterDerivative.getSize<0>()
        );
        MatrixMapBlock projectedView(
            projectedMap,
            0, 0,
            _projectedParameterDerivative.getSize<1>(), 
            _projectedParameterDerivative.getSize<0>()
        );
        // END TODO
        morphologyView = projectedView * _morphologyProjection->computeProjectedParameterJacobian();
    }
}

void multifit::ComponentModelProjection::_computeWcsParameterDerivative(
    ndarray::Array<Pixel,2,1> const & matrix
) {
    if (hasTranslationDerivative()) {
        _ensureTranslationDerivative();
        MatrixMap astrometryMap(
            matrix.getData(),
            matrix.getStride<0>(),
            getModel()->getAstrometry()->getParameterSize()
        );                      
        MatrixMapBlock astrometryView(
            astrometryMap, 
            0, 0, 
            matrix.getSize<1>(), matrix.getSize<0>()
        );
        astrometryView = getTranslationMatrixView() * _transform->dTransform(
            getModel()->getAstrometry()->computePosition()
        );
    }
    if (hasProjectedParameterDerivative()) {
        _ensureProjectedParameterDerivative();
        // TODO: Move this into an inline function when possible.
        int offset = getModel()->getAstrometry()->getParameterSize();
        MatrixMap morphologyMap(
            matrix.getData(),
            matrix.getStride<0>(), 
            getModel()->getMorphology()->getNonlinearParameterSize()
        );
        MatrixMapBlock morphologyView(
            morphologyMap, 
            0, offset, 
            matrix.getSize<1>(), morphologyMap.cols()
        );
        // END TODO
        // TODO: Move this into an inline function when possible.
        MatrixMap projectedMap(
            _projectedParameterDerivative.getData(),
            _projectedParameterDerivative.getStride<0>(), 
            _projectedParameterDerivative.getSize<0>()
        );
        MatrixMapBlock projectedView(
            projectedMap,
            0, 0,
            _projectedParameterDerivative.getSize<1>(), 
            _projectedParameterDerivative.getSize<0>()
        );
        // END TODO
	morphologyView += projectedView *
            _morphologyProjection->computeTransformParameterJacobian();
    }
}

/**
 * Ensure's that _translationDerivative is up to date
 *
 * If _translationDerivative array has not been allocated, do so first.
 * If it is not up to date, call _computeTranslationParameterDerivative
 */
void multifit::ComponentModelProjection::_ensureTranslationDerivative() {
    if (_translationDerivative.empty()) {
        ndarray::shallow(_translationDerivative) = ndarray::allocate<Allocator>(
            ndarray::makeVector(2, getFootprint()->getNpix())
        );
        _validProducts &= (~TRANSLATION_DERIVATIVE);
    }
    if (!(_validProducts & TRANSLATION_DERIVATIVE)) {
        _computeTranslationDerivative(_translationDerivative);
    }
}

/**
 * Ensure that _projectedParameterDerivative is up to date
 *
 * If _projectedParameterDerivative array has not been allocated, do so first.
 * If it is not up to date, call _computeProjectedParameterDerivative
 */
void multifit::ComponentModelProjection::_ensureProjectedParameterDerivative() {
    if (_projectedParameterDerivative.empty()) {
        ndarray::shallow(_projectedParameterDerivative) =  ndarray::allocate<Allocator>(
                ndarray::makeVector(2, getFootprint()->getNpix())
        );
        _validProducts &= (~PROJECTED_PARAMETER_DERIVATIVE);
    }
    if (!(_validProducts & PROJECTED_PARAMETER_DERIVATIVE)) {
        _computeProjectedParameterDerivative(_projectedParameterDerivative);
    }
}
