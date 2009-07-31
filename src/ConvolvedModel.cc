#include "lsst/meas/multifit/ConvolvedModel.h"

namespace multifit = lsst::meas::multifit;

void multifit::ConvolvedModel::init() {
    _psfImage.resize(getImageSize());
    _dPsfdX.resize(getImageSize());
    _dPsfDy.resize(getImageSize());
    _imageDirty = true;
    _linearDirty = true; 
    _nonlinearDirty = true;
    _psfDirty = true;
    _transformDirty = true;
}

void multifit::ConvolvedModel::updatePsfProducts() {
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

void multifit::ConvolvedModel::updateParametrizedImage() {
    _unconvolved->updateParametrizedImage();
    ImageVector convolvedRef = getImageView(_parametrizedImage,
            getImageHeight(), getImageWidth()
    );
    ImageVector unconvolvedRef = getImageView(_unconvolved->_parametrizedImage,
            _unconvolved->getImageHeight(), _unconvolved->getImageWidth()
    );
    _psf->getKernel()->convolve(unconvolvedRef, convolvedRef, oversample);
}


void multifit::ConvolvedModel::updateConstantImage() {
    if(!hasConstantImage())
        return;

    _unconvolved->updateConstantImage();
    ImageVector convolvedRef = getImageView(_constantImage,
            getImageHeight(), getImageWidth()
    );
    ImageVector unconvolvedRef = getImageView(_unconvolved->_constantImage,
            _unconvolved->getImageHeight(), _unconvolved->getImageWidth()
    );
    _psf->getKernel()->convolve(unconvolvedRef, convolvedRef, oversample);
}

void multifit::ConvolvedModel::updateLinearMatrix() {
    _unconvolved->updateLinearMatrix();
    for(int i = 0; i < getLinearSize(); ++i) {
        ImageVector convolvedRef = getImageView(_linearMatrix.row(i),
               getImageHeight(), getImageWidth()
        );
        ImageVector unconvolvedRef = getImageView(
                _unconvolved->_linearMatrix.row(i),
                _unconvolved->getImageHeight(), _unconvolved->getImageWidth()
        );
        _psf->getKernel()->convolve(unconvolvedRef, convolvedRef, oversample);
    }
}

void multifit::ConvolvedModel::updateNonlinearMatrix() {
    _unconvolved->updateNonlinearMatrix();
    for(int i = 0; i < getNonLinearSize(); ++i) {
        ImageVector convolvedRef = getImageView(_nonlinearMatrix.row(i),
               getImageHeight(), getImageWidth()
        );
        ImageVector unconvolvedRef = getImageView(
                _unconvolved->_nonlinearMatrix.row(i),
                _unconvolved->getImageHeight(), _unconvolved->getImageWidth()
        );
        _psf->getKernel()->convolve(unconvolvedRef, convolvedRef, oversample);
    }
}

void multifit::ConvolvedModel::updateTransformMatrix() {
    _unconvolved->updateTransformMatrix();
    //TODO:: implement me
}

void multifit::ConvolvedModel::updatePsfMatrix() {
    //TODO:: implement me
}

