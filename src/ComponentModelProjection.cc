#include "lsst/meas/multifit/ComponentModelProjection.h"
#include "lsst/meas/multifit/matrices.h"

namespace multifit = lsst::meas::multifit;

multifit::ComponentModelProjection::ComponentModelProjection(
    ComponentModel::ConstPtr const & model,
    KernelConstPtr const & kernel,
    WcsConstPtr const & wcs,
    FootprintConstPtr const & footprint,
    double photFactor
) : ModelProjection(model, wcs, footprint, photFactor),
    _validProducts(0), 
    //TODO: need wcs linearize api
    //_transform(wcs->linearize(model->getAstrometry()->apply())),
    _morphologyProjection(), 
    _translationDerivative(), 
    _projectedParameterDerivative()
{
    _morphologyProjection = model->getMorphology()->makeProjection(
        lsst::afw::geom::ExtentI::makeXY(
            kernel->getWidth(), kernel->getHeight()
        ), 
        _transform
    );
}

void multifit::ComponentModelProjection::_computeNonlinearParameterDerivative(
    ndarray::Array<Pixel,2,2> const & matrix
) {
    _ensureTranslationDerivative();
    getAstrometryParameterMatrixView(matrix) = getTranslationMatrixView() * 
        _transform->matrix().linear() * getModel()->getAstrometry()->differentiate();
    _ensureProjectedParameterDerivative();
    getMorphologyParameterMatrixView(matrix) = getProjectedParameterMatrixView() *
        _morphologyProjection->computeProjectedParameterJacobian();
}

void multifit::ComponentModelProjection::_computeWcsParameterDerivative(
    ndarray::Array<Pixel,2,2> const & matrix
) {
    _ensureTranslationDerivative();
    getAstrometryParameterMatrixView(matrix) = getTranslationMatrixView() *
        _transform->dTransform(getModel()->getAstrometry()->apply());
    _ensureProjectedParameterDerivative();
    getMorphologyParameterMatrixView(matrix) += getProjectedParameterMatrixView() *
        _morphologyProjection->computeTransformParameterJacobian();    
}

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
