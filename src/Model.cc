#include "lsst/meas/multifit/Model.h"

namespace multifit = lsst::meas::multifit;

void init(int const & nonlinearSize,
        int const & linearsize,
        int const & psfBasisSize) {
    int nPix = getImageSize();
    if(nPix != 0) {
        _parameterizedImage(nPix);
        _constantImage(nPix);
        
        if(linearSize != 0) {
            _linearMatrix.resize(linearSize, nPix);
        }
        if(nonlinearSize != 0) {
            _nonlinearMatrix.resize(nonlinearSize, nPix);    
        }
        if(psfBasisSize != 0)
            _psfMatrix.resize(psfBasisSize, nPix);
        
        _transformMatrix.resize(TRANSFORM_SIZE, nPix);
    }   
}

Eigen::VectorXd const & multifit::Model::computeParameterizedImage() {
    updateParameterizedImage();
    return _parameterizedImage;
}

Eigen::VectorXd const & multifit::Model::computeConstantImage() {
    updateConstantImage();
    return _constantImage;
}

Eigen::MatrixXd const & multifit::Model::computeLinearMatrix()  {
    updateLinearMatrix();
    return _linearMatrix;
}

Eigen::MatrixXd const & multifit::Model::computeNonlinearMatrix() {
    updateNonlinearMatrix();
    return _nonlinearMatrix;
}

Eigen::MatrixXd const & multifit::Model::computeTransformMatrix() {
    updateTransformMatrix();
    return _transformMatrix;
}

Eigen::MatrixXd const & multifit::Model::computePsfMatrix() {
    updatePsfMatrix();
    return _psfMatrix();
}
