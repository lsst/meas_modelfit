#include "lsst/meas/multifit/PointSourceModel.h"
#include "lsst/afw/image/Image.h"

namespace multifit = lsst::meas::multifit;

void multifit::PointSourceModel::evalParametrizedImage(
        ImageVector const & modelImage
) const {
    if (!_psf) {       
        throw exception("Null psf in PointSourceModel");
    }
    ImageVector::Index psfShape = modelImage.shape();
    ImageVector::Index psfStride(psfShape[1], 1);
    double * psfData = new double[psfShape.products()];
    ImageVector psfImage(psfData, psfShape, psfStride);
    _psf->getKernel()->getImage(psfImage, _center);

    ImageVector::iterator imageIter = modelImage.begin();
    ImageVector::iterator const imageEnd = modelImage.end();

    ImageVector::iterator psfIter = psfImage->begin();

    for( ; imageIter != imageEnd; ++imageIter, ++psfIter) {
        *imageIter += *psfImage * _magnitude;
    }

    delete psfData;
}

void multifit::PointSourceModel::evalLinearDerivative(
        DerivativeMatrix const & linearDerivative
) const {
    if (!_psf) {       
        throw exception("Null psf in PointSourceModel");
    }
    
    _psf->getKernel()->getImage(linearDerivative[0], _center);
}

void multifit::PointSourceModel::evalNonlinearDerivative(
        DerivativeMatrix const & nonlinearDerivative
) const {
    ImageVector::Index imageShape = nonlinearDerivatives[0].shape();
    int nPix = imageShape.products();

    Eigen::VectorXd dPsfDx(nPix), dPsfDy(nPix);
    ndarray::ArrafRef<double, 2, 1> dPsfDxRef(dPsfDx.data(), imageShape);
    ndarray::ArrafRef<double, 2, 1> dPsfDyRef(dPsfDy.data(), imageShape);
    _psf->getKernel()->getCentroidDerivative(dPsfDxRef, dPsfDyRef);

    Eigen::Map<VectorXd> dxView = extractEigenView(nonlinearDerivative[0]);
    Eigen::Map<VectorXd> dyView = extractEigenView(nonlinearDerivative[1]);
    
    dxView = dPsfDx * _transform(0, 0) + dPsfDy * _transform(0, 1);
    dyView = dPsfDx * _transform(1, 0) + dPsfDy * _transform(1, 1);        
}

void multifit::PointSourceModel::evalTransformDerivative(
        DerivativeMatrix const & transformDerivative
) const {
    /TODO: implement me            
}

void multifit::PointSourceModel::evalPsfDerivative(
        DerivativeMatrix const & psfDerivative     
) const {
        
}
