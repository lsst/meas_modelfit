#include <lsst/meas/multifit/ModelProjection.h>
#include <lsst/meas/multifit/Model.h>
#include <lsst/meas/multifit/matrices.h>

namespace projections  = lsst::meas::multifit::projections;

ndarray::Array<Pixel const,1,1> projections::ModelProjection::computeModelImage() {
    if(!(_activeProducts & ModelProjection::MODEL_IMAGE)) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Product MODEL_IMAGE is not enabled."
        );
    }
    if (_modelImage.empty()) {
        ndarray::shallow(_modelImage) = ndarray::allocate<Allocator>(
            ndarray::makeVector(_footprint->getNpix())
        );
        _validProducts &= (~MODEL_IMAGE);
    }
    if (!(_validProducts & MODEL_IMAGE)) {
        _computeModelImage(_modelImage);
        _validProducts |= MODEL_IMAGE;
    }
    return _modelImage;
}

ndarray::Array<Pixel const,2,2> projections::ModelProjection::computeLinearParameterDerivative() {
    if(!(_activeProducts & ModelProjection::LINEAR_PARAMETER_DERIVATIVE)) {     
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Product MODEL_IMAGE is not enabled."
        );
    }
    if (_linearParameterDerivative.empty()) {
        ndarray::shallow(_linearParameterDerivative) = ndarray::allocate<Allocator>(
            ndarray::makeVector(getLinearParameterSize(),_footprint->getNpix())
        );
        _validProducts &= (~LINEAR_PARAMETER_DERIVATIVE);
    }
    if (!(_validProducts & LINEAR_PARAMETER_DERIVATIVE)) {
        _computeLinearParameterDerivative(_linearParameterDerivative);
        _validProducts |= LINEAR_PARAMETER_DERIVATIVE;
    }
    return _linearParameterDerivative;
}

ndarray::Array<Pixel const,2,2> projections::ModelProjection::computeNonlinearParameterDerivative() {
    if (!(_activeProducts & NONLINEAR_PARAMETER_DERIVATIVE)) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Product NONLINEAR_PARAMETER_DERIVATIVE is not enabled."
        );
    }
    if (_nonlinearParameterDerivative.empty()) {
        ndarray::shallow(_nonlinearParameterDerivative) = ndarray::allocate<Allocator>(
            ndarray::makeVector(
                getNonlinearParameterSize(),
                _footprint->getNpix()
            )
        );
        _validProducts &= (~NONLINEAR_PARAMETER_DERIVATIVE);
    }
    if (!(_validProducts & NONLINEAR_PARAMETER_DERIVATIVE)) {
        _computeNonlinearParameterDerivative(_nonlinearParameterDerivative);
        _validProducts |= NONLINEAR_PARAMETER_DERIVATIVE;
    }
    return _nonlinearParameterDerivative;
}

ndarray::Array<Pixel const,2,2> projections::ModelProjection::computeWcsParameterDerivative() {
    if (!(_activeProducts & WCS_PARAMETER_DERIVATIVE)) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Product WCS_PARAMETER_DERIVATIVE is not enabled."
        );
    }
    if (_wcsParameterDerivative.empty()) {
        ndarray::shallow(_wcsParameterDerivative) = ndarray::allocate<Allocator>(
            ndarray::makeVector(getWcsParameterSize(),_footprint->getNpix())
        );
        _validProducts &= (~WCS_PARAMETER_DERIVATIVE);
    }
    if (!(_validProducts & WCS_PARAMETER_DERIVATIVE)) {
        _computeWcsParameterDerivative(_wcsParameterDerivative);
        _validProducts |= WCS_PARAMETER_DERIVATIVE;
    }
    return _wcsParameterDerivative;
}

ndarray::Array<Pixel const,2,2> projections::ModelProjection::computePsfParameterDerivative() {
    if(!(_activeProducts & ModelProjection::PSF_PARAMETER_DERIVATIVE)) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Product PSF_PARAMETER_DERIVATIVE is not enabled."
        );    
    }
    if (_psfParameterDerivative.empty()) {
        ndarray::shallow(_psfParameterDerivative) = ndarray::allocate<Allocator>(
            ndarray::makeVector(getPsfParameterSize(),_footprint->getNpix())
        );
        _validProducts &= (~PSF_PARAMETER_DERIVATIVE);
    }
    if (!(_validProducts & PSF_PARAMETER_DERIVATIVE)) {
        _computePsfParameterDerivative(_psfParameterDerivative);
        _validProducts |= PSF_PARAMETER_DERIVATIVE;
    }
    return _psfParameterDerivative;
}

void projections::ModelProjection::setModelImageBuffer(
    ndarray::Array<Pixel,1,1> const & buffer
) {
    if (buffer.getSize<0>() != _footprint->getNpix()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptioins::LengthErrorException,
            "Model image buffer's size must match Footprint size."
        );
    }
    ndarray::shallow(_modelImage) = buffer;
    _validProducts &= (~MODEL_IMAGE);
}

void projections::ModelProjection::setLinearParameterDerivativeBuffer(
    ndarray::Array<Pixel,2,2> const & buffer
) {
    if (buffer.getSize<1>() != _footprint->getNpix()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LengthErrorException,
            "Linear parameter derivative buffer's image dimension "
            "must match Footprint size."
        );
    }
    if (buffer.getSize<0>() != _model->getLinearParameterSize()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LengthErrorException,
            "Linear parameter derivative buffer's parameter dimension "
            "must match linear parameter size."
        );
    }
    ndarray::shallow(_linearParameterDerivative) = buffer;
    _validProducts &= (~LINEAR_PARAMETER_DERIVATIVE);
}

void projections::ModelProjection::setNonlinearParameterDerivativeBuffer(
    ndarray::Array<Pixel,2,2> const & buffer
) {
    if (buffer.getSize<1>() != _footprint->getNpix()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LengthErrorException,
            "Nonlinear parameter derivative buffer's image dimension "
            "must match Footprint size."
        );
    }
    if (buffer.getSize<0>() != _model->getNonlinearParameterSize()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LengthErrorException,
            "Nonlinear parameter derivative buffer's parameter dimension "
            "must match nonlinear parameter size."
        );
    }
    ndarray::shallow(_nonlinearParameterDerivative) = buffer;
    _validProducts &= (~NONLINEAR_PARAMETER_DERIVATIVE);
}

void projections::ModelProjection::setWcsParameterDerivativeBuffer(
    ndarray::Array<Pixel,2,2> const & buffer
) {
    if (buffer.getSize<1>() != _footprint->getNpix()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LengthErrorException,
            "WCS parameter derivative buffer's image dimension "
            "must match Footprint size."
        );
    }
    if (buffer.getSize<0>() != getWcsParameterSize()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LengthErrorException,
            "WCS parameter derivative buffer's parameter dimension "
            "must match wcs parameter size."
        );
    }
    ndarray::shallow(_wcsParameterDerivative) = buffer;
    _validProducts &= (~WCS_PARAMETER_DERIVATIVE);
}

void projections::ModelProjection::setPsfParameterDerivativeBuffer(
    ndarray::Array<Pixel,2,2> const & buffer
) {
    if (buffer.getSize<1>() != _footprint->getNpix()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LengthErrorException,
            "PSF parameter derivative buffer's image dimension "
            "must match Footprint size."
        );
    }
    if (buffer.getSize<0>() != getPsfParameterSize()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LengthErrorException,
            "PSF parameter derivative buffer's parameter dimension "
            "must match psf parameter size."
        );
    }
    ndarray::shallow(_psfParameterDerivative) = buffer;
    _validProducts &= (~PSF_PARAMETER_DERIVATIVE);
}

void projections::ModelProjection::_computeModelImage(
    ndarray::Array<Pixel,1,1> const & vector
) {
    getCompressedVectorView(vector) = 
        getCompressedMatrixView(computeLinearParameterDerivative()) * 
        getLinearParameters();
}

void projections::ModelProjection::_handleLinearParameterChange() {
    _validProducts &= (~MODEL_IMAGE);
    _validProducts &= (~NONLINEAR_PARAMETER_DERIVATIVE);
    _validProducts &= (~WCS_PARAMETER_DERIVATIVE);
    _validProducts &= (~PSF_PARAMETER_DERIVATIVE);
}

void projections::ModelProjection::_handleNonlinearParameterChange() {
    _validProducts &= (~MODEL_IMAGE);
    _validProducts &= (~LINEAR_PARAMETER_DERIVATIVE);
    _validProducts &= (~NONLINEAR_PARAMETER_DERIVATIVE);
    _validProducts &= (~WCS_PARAMETER_DERIVATIVE);
    _validProducts &= (~PSF_PARAMETER_DERIVATIVE);
}
