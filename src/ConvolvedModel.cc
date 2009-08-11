#include "lsst/meas/multifit/ConvolvedModel.h"

namespace multifit = lsst::meas::multifit;

void multifit::ConvolvedModel::init() {
    _constantDirty = true;
    _parametrizedDirty = true;
    _linearDirty = true; 
    _nonlinearDirty = true;
    _psfDirty = true;
    _transformDirty = true;
}

void multifit::ConvolvedModel::updateParametrizedImage() {
    if(!_parametrizedDirty)
        return;

    VectorXd const & unconvolvedImage = _unconvolved->computeParametrizedImage();

    ImageVector outputRef = getImageView(_parametrizedImage,
            getImageHeight(), getImageWidth());
    ImageVector inputRef = getImageView(unconvolvedImage,
            _unconvolved->getImageHeight(), _unconvolved->getImageWidth());
    
    _psf->getKernel()->convolve(inputRef, outputRef, _oversample);
    _parametrizedDirty = false;
}


void multifit::ConvolvedModel::updateConstantImage() {
    if(!hasConstantImage() && !_constantDirty)
        return;

    VectorXd const & unconvolvedImage = _unconvolved->computeConstantImage();

    ImageVector outputRef = getImageView(_constantImage,
            getImageHeight(), getImageWidth());
    ImageVector inputRef = getImageView(unconvolvedImage,
            _unconvolved->getImageHeight(), _unconvolved->getImageWidth());

    _psf->getKernel()->convolve(inputRef, ouputRef, _oversample);
    _constantDirty = false;
}

void multifit::ConvolvedModel::updateLinearMatrix() {
    if(!_linearDirty)
        return;

    MatrixXd const & unconvolvedLinear = _unconvolved->computeLinearMatrix();

    int height = getImageHeight();
    int width = getImageWidth();
    int unconvolvedHeight = _unconvolved->getImageHeight();
    int unconvolvedWidth = _unconvolved->getImageWidth();

    for(int i = 0; i < getLinearSize(); ++i) {
        ImageVector outputRef = getImageView(_linearMatrix, i,
                height, width);
        ImageVector inputRef = getImageView(unconvolvedLinear, i,
                unconvolvedHeight, unconvolvedWidth);
        _psf->getKernel()->convolve(inputRef, outputRef, _oversample);
    }
    _linearDirty = false;
}

void multifit::ConvolvedModel::updateNonlinearMatrix() {
    if(!_nonlinearDirty)
        return;

    MatrixXd const & unconvolvedNonlinear = _unconvolved->computeNonlinearMatrix();

    int height = getImageHeight();
    int width = getImageWidth();
    int unconvolvedHeight = _unconvolved->getImageHeight();
    int unconvolvedWidth = _unconvolved->getImageWidth();
    
    for(int i = 0; i < getNonLinearSize(); ++i) {
        ImageVector ouputRef = getImageView(_nonlinearMatrix, i,
                height, width);
        ImageVector inputRef = getImageView(unconvolvedNonlinear, i,
                unconvolvedHeight, unconvolvedWidth);
        _psf->getKernel()->convolve(inputRef, outputRef, _oversample);
    }

    _nonlinearDirty = false;
}

void multifit::ConvolvedModel::updateTransformMatrix() {
    if(!_transformDirty)
        return;

    MatrixXd const & unconvolvedTranform = _unconvolved->computeTransformMatrix();
    
    int height = getImageHeight();
    int width = getImageWidth();
    int unconvolvedHeight = _unconvolved->getImageHeight();
    int unconvolvedWidth = _unconvolved->getImageWidth();

    for(int i = 0; i < getPsfBasisSize(); ++i) {
        ImageVector outputRef = getImageView(_transformMatrix, i,
                height, width);
        ImageVector inputRef = getImageView(_unconvolvedTransform, i,
                unconvolvedHeight, unconvolvedWidth);
        _psf->getKernel()->convolve(inputRef, outputRef, _oversample);
    }

    _transformDirty = false;
}

void multifit::ConvolvedModel::updatePsfMatrix() {
    if(_psfDirty)
        return;

    int unconvolvedHeight = _unconvolved->getImageHeight();
    int unconvolvedWidth = _unconvolved->getImageWidth();

    VectorXd unconvolvedImage;
    if(_parametrizedDirty)
        unconvolvedImage << _unconvolved->computeParmetrizedImage();
    else unconvolvedImage << _unconvolved->_parametrizedImage;
    if (hasConstantImage()){
        if(_constantDirty)
            unconvolvedImage += _unconvoled->computeConstantImage();
        else unconvolvedImage += _unconvolved->_constantImage;
    }

    ImageVector inputRef = getImageView(unconvolvedImage,
            unconvolvedHeight, unconvolvedWidth
    );

    int height = getImageHeight();
    int width = getImageWidth();

    for(int i = 0; i < getPsfBasisSize(); ++i) {
        ImageVector outputRef = getImageView(_psfMatrix, i,
            height, width);
        _psf->getBasisKernel(i)->convolve(inputRef, outputRef, _oversample);
    }

    _psfDirty = false;
}

