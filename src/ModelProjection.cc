#include "lsst/meas/multifit/ModelProjection.h"
#include "lsst/meas/multifit/matrices.h"
#include "lsst/pex/exceptions/Runtime.h"

namespace multifit = lsst::meas::multifit;

multifit::ModelProjection::ModelProjection(
    Model::ConstPtr const & model,
    WcsConstPtr const & wcs,
    FootprintConstPtr const & footprint
) : _validProducts(0),
    _model(model),
    _footprint(footprint), _wcs(wcs),
    _modelImage(), 
    _linearParameterDerivative(), 
    _nonlinearParameterDerivative(),
    _wcsParameterDerivative(), 
    _psfParameterDerivative()
{
    if (_footprint->getNpix() <= 0) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Cannot create model projection with empty footprint."
        );
    }
}

ndarray::Array<multifit::Pixel const,1,1> multifit::ModelProjection::computeModelImage() {
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

ndarray::Array<multifit::Pixel const,2,1> multifit::ModelProjection::computeLinearParameterDerivative() {
    if (!hasLinearParameterDerivative()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Cannot compute linear parameter derivative."
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

ndarray::Array<multifit::Pixel const,2,1> multifit::ModelProjection::computeNonlinearParameterDerivative() {
    if (!hasNonlinearParameterDerivative()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Cannot compute nonlinear parameter derivative."
        );
    }
    if (_nonlinearParameterDerivative.empty()) {
        ndarray::shallow(_nonlinearParameterDerivative) = ndarray::allocate<Allocator>(
            ndarray::makeVector(getNonlinearParameterSize(),_footprint->getNpix())
        );
        _validProducts &= (~NONLINEAR_PARAMETER_DERIVATIVE);
    }
    if (!(_validProducts & NONLINEAR_PARAMETER_DERIVATIVE)) {
        _computeNonlinearParameterDerivative(_nonlinearParameterDerivative);
        _validProducts |= NONLINEAR_PARAMETER_DERIVATIVE;
    }
    return _nonlinearParameterDerivative;
}

ndarray::Array<multifit::Pixel const,2,1> multifit::ModelProjection::computeWcsParameterDerivative() {
    if (!hasWcsParameterDerivative()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Cannot compute wcs parameter derivative."
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

ndarray::Array<multifit::Pixel const,2,1> multifit::ModelProjection::computePsfParameterDerivative() {
    if (!hasPsfParameterDerivative()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Cannot compute psf parameter derivative."
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

void multifit::ModelProjection::setModelImageBuffer(ndarray::Array<Pixel,1,1> const & buffer) {
    if (buffer.getSize<0>() != _footprint->getNpix()) {
        throw std::invalid_argument("Model image buffer's size must match Footprint size.");
    }
    ndarray::shallow(_modelImage) = buffer;
    _validProducts &= (~MODEL_IMAGE);
}

void multifit::ModelProjection::setLinearParameterDerivativeBuffer(ndarray::Array<Pixel,2,1> const & buffer) {
    if (!hasLinearParameterDerivative()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Projection has no linear parameter derivative."
        );
    }
    if (buffer.getSize<1>() != _footprint->getNpix()) {
        throw std::invalid_argument("Linear parameter derivative buffer's image dimension "
                                    "must match Footprint size.");
    }
    if (buffer.getSize<0>() != _model->getLinearParameterSize()) {
        throw std::invalid_argument("Linear parameter derivative buffer's parameter dimension "
                                    "must match linear parameter size.");
    }
    ndarray::shallow(_linearParameterDerivative) = buffer;
    _validProducts &= (~LINEAR_PARAMETER_DERIVATIVE);
}

void multifit::ModelProjection::setNonlinearParameterDerivativeBuffer(
    ndarray::Array<Pixel,2,1> const & buffer
) {
    if (!hasNonlinearParameterDerivative()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Projection has no nonlinear parameter derivative."
        );
    }
    if (buffer.getSize<1>() != _footprint->getNpix()) {
        throw std::invalid_argument("Nonlinear parameter derivative buffer's image dimension "
                                    "must match Footprint size.");
    }
    if (buffer.getSize<0>() != _model->getNonlinearParameterSize()) {
        throw std::invalid_argument("Nonlinear parameter derivative buffer's parameter dimension "
                                    "must match nonlinear parameter size.");
    }
    ndarray::shallow(_nonlinearParameterDerivative) = buffer;
    _validProducts &= (~NONLINEAR_PARAMETER_DERIVATIVE);
}

void multifit::ModelProjection::setWcsParameterDerivativeBuffer(ndarray::Array<Pixel,2,1> const & buffer) {
    if (!hasWcsParameterDerivative()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Projection has no wcs parameter derivative."
        );
    }
    if (buffer.getSize<1>() != _footprint->getNpix()) {
        throw std::invalid_argument("WCS parameter derivative buffer's image dimension "
                                    "must match Footprint size.");
    }
    if (buffer.getSize<0>() != getWcsParameterSize()) {
        throw std::invalid_argument("WCS parameter derivative buffer's parameter dimension "
                                    "must match WCS parameter size.");
    }
    ndarray::shallow(_wcsParameterDerivative) = buffer;
    _validProducts &= (~WCS_PARAMETER_DERIVATIVE);
}

void multifit::ModelProjection::setPsfParameterDerivativeBuffer(ndarray::Array<Pixel,2,1> const & buffer) {
    if (!hasPsfParameterDerivative()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Projection has no psf parameter derivative."
        );
    }
    if (buffer.getSize<1>() != _footprint->getNpix()) {
        throw std::invalid_argument("PSF parameter derivative buffer's image dimension "
                                    "must match Footprint size.");
    }
    if (buffer.getSize<0>() != getPsfParameterSize()) {
        throw std::invalid_argument("PSF parameter derivative buffer's parameter dimension "
                                    "must match PSF parameter size.");
    }
    ndarray::shallow(_psfParameterDerivative) = buffer;
    _validProducts &= (~PSF_PARAMETER_DERIVATIVE);
}

void multifit::ModelProjection::_computeModelImage(ndarray::Array<Pixel,1,1> const & vector) {
    ndarray::Array<Pixel const,2,1> array(computeLinearParameterDerivative());
    MatrixMap map(array.getData(), array.getStride<0>(), array.getSize<0>());
    MatrixMapBlock block(map, 0, 0, array.getSize<1>(), array.getSize<0>());
    getVectorView(vector) = block * _model->getLinearParameterVector();
}

void multifit::ModelProjection::_handleLinearParameterChange() {
    _validProducts &= (~MODEL_IMAGE);
    _validProducts &= (~NONLINEAR_PARAMETER_DERIVATIVE);
    _validProducts &= (~WCS_PARAMETER_DERIVATIVE);
    _validProducts &= (~PSF_PARAMETER_DERIVATIVE);
}

void multifit::ModelProjection::_handleNonlinearParameterChange() {
    _validProducts &= (~MODEL_IMAGE);
    _validProducts &= (~LINEAR_PARAMETER_DERIVATIVE);
    _validProducts &= (~NONLINEAR_PARAMETER_DERIVATIVE);
    _validProducts &= (~WCS_PARAMETER_DERIVATIVE);
    _validProducts &= (~PSF_PARAMETER_DERIVATIVE);
}

