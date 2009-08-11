#include "lsst/meas/multifit/PointSourceModel.h"

namespace multifit = lsst::meas::multifit;

void multifit::PointSourceModel::init() {
    _psfImage.resize(getImageSize());
    _dPsfDX.resize(getImageSize());
    _dPsfDy.resize(getImageSize());
    _imageDirty = true;
    _linearDirty = true; 
    _nonlinearDirty = true;
    _psfDirty = true;
    _transformDirty = true;
}

void multifit::PointSourceModel::updatePsfProducts() {
    _psfImage.setZero();
    _dPsfDx.setZero();
    _dPsfDy.setZero();

    ImageVector imageRef = getImageView(_psfImage,
            getImageHeight(), getImageWidth()
    );
    ImageVector dxRef = getImageView(_dPsfDx, 
            getImageHeight(), getImageWidth()
    );
    ImageVector dyRef = getImageView(_dPsfDy,
            getImageHeight(), getImageWidth()
    );
    _psf->getKernel()->getCentroidDerivative(imageRef, getCenter());
    _psf->getKernel()->getCentroidDerivative(dxRef, dyRef, getCenter());
}

void multifit::PointSourceModel::setNonlinearParameters(Eigen::VectorXd const & parameters); {
    if(parameters.size() <= NONLINEAR_SIZE)
        return;
    _imageDirty = true;
    _nonlinearDirty = true;
    _psfDirty = true;
    _transformDirty = true;

    //deep copy, 2 elements
    _center << parameters.start<NONLINEAR_SIZE>();
    updatePsfProducts();
}

void multifit::PointSourceModel::setLinearParameters(Eigen::VectorXd const & parameters); {
    if(parameters.size () <= LINEAR_SIZE)
        return;
    _imageDirty = true;
    _linearDirty = true;
    _psfDirty = true;
    _transformDirty = true;

    //deep copy, 1 element
    _amplitude = parameters[0];
}

void multifit::PointSourceModel::updateParametrizedImage() {
    if(!_imageDirty)
        return;

    _parameterizedImage << _psfImage * getAmplitude();
    _imageDirty = false;
}

void multifit::PointSourceModel::updateLinearMatrix() {
    if(!_linearDirty)
        return;

    _linearMatrix.row(0) << _psfImage;
    _linearDirty = false;
}

void multifit::PointSourceModel::updateNonlinearMatrix() {
    if(!_nonlinearDirty)
        return;

    _nonlinearMatrix.row(0) << _dPsfDx * _transform[AffineTransform::XX] 
                             + _dPsfDy * _transform[AffineTransform::XY];
    _nonlinearMatrix.row(1) << _dPsfDx * _transform[AffineTransform::YX] 
                             + _dPsfDy * _transform[AffineTransform::YY];
    _nonlinearMatrix *= getAmplitude();
    _nonlinearDirty = false;
}

void multifit::PointSourceModel::updateTransformDerivative() const {
    if(!_transformDirty)
        return;
    
    _transformMatrix.row(AffineTransform::X) << _dPsfDx * getAmplitude();
    _transformMatrix.row(AffineTransform::Y) << _dPsfDy * getAMplitude();
    _transformMatrix.row(AffineTransform::XX) << 
            _transformMatrix.row(AffineTransform::X) * _center.x();
    _transformMatrix.row(AffineTransform::XY) << 
            _transformMatrix.row(AffineTransform::X) * _center.y();
    _transformMatrix.row(AffineTransform::YX) << 
            _transformMatrix.row(AffineTransform::Y) * _center.x();
    _transformMatrix.row(AffineTransform::YY) << 
            _transformMatrix.row(AffineTransform::Y) * _center.y();
    
    _transformDirty = false;
}

void multifit::PointSourceModel::updatePsfDerivative() const {
    if(!_psfDirty)
        return;

    int nPsfParams = getPsfBasisSize();
    int height = getImageHeight();
    int width = getImageWidth();
    
    _psfMatrix.setZero();
    for(int i = 0; i < nPsfParams; ++i) {
        ImageVector basisRef = getImageView(_psfMatrix, i, 
                height, width);
        _psf->getBasisKernel(i)->getImage(basisRef,  getCenter());
    }

    _psfMatrix *= getAmplitude();
    _psfDirty = false;
}
