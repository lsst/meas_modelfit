#include "lsst/meas/multifit/PointSourceModel.h"
#include "lsst/afw/image/Image.h"

namespace multifit = lsst::meas::multifit;

void multifit::PointSourceModel::evalParametrizedImage(
        ImageVector const & modelImage
) const {
    ImageVector::Index psfShape = modelImage.shape();
    int nPix = psfShape.products();
    Eigen::VectorXd psfImage(nPix);
    ndarray::ArrayRef<double, 2, 1> psfRef(psfImage.data(), psfShape);
    _psf->getKernel()->getImage(psfRef, _center);

    Eigen::Map<VectorXd> imageView = extractEigenView(modelImage);
    imageView << psfImage * _magnitude
}

void multifit::PointSourceModel::evalLinearDerivative(
        DerivativeMatrix const & linearDerivative
) const {
    assert(linearDerivative.size() == getNumLinearParameters());

    _psf->getKernel()->getImage(linearDerivative[0], _center);
}

void multifit::PointSourceModel::evalNonlinearDerivative(
        DerivativeMatrix const & nonlinearDerivative
) const {
    assert(nonlinearDerivative.size() == getNumNonlinearParameters());

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

    dxView *= _magnitude;
    dyView *= _magnitude;
}

void multifit::PointSourceModel::evalTransformDerivative(
        DerivativeMatrix const & transformDerivative
) const {
    assert(transformDerivative.size() == getNumTransformParameters());

    _psf->getKernel()->getCentroidDerivative(
            transformDerivative[4],
            transformDerivative[5],
            _center
    );

    Eigen::Map<VectorXd> dTxView = extractEigenView(transformDerivative[4]);
    Eigen::Map<VectorXd> dTyView = extractEigenView(transformDerivative[5]);
    
    Eigen::Map<VectorXd> dTxxView = extractEigenView(transformDerivative[0]);
    Eigen::Map<VectorXd> dTxyView = extractEigenView(transformDerivative[1]);
    Eigen::Map<VectorXd> dTyxView = extractEigenView(transformDerivative[2]);
    Eigen::Map<VectorXd> dTyyView = extractEigenView(transformDerivative[3]);


    dTxView *= _magnitude;
    dTyView *= _magnitude;
    dTxxView << dTxView * _center.x();
    dTxyView << dTxView * _center.y();
    dTyxView << dTyView * _center.x();
    dTyyView << dTyView * _center.y();
}

void multifit::PointSourceModel::evalPsfDerivative(
        DerivativeMatrix const & psfDerivative     
) const {
    int nPsfParams = getNumPsfParameters();
    assert(psfDerivative.size() == nPsfParams);

    Eigen::VectorXd coefficients = _psf->getCoefficients();

    for(int i = 0; i < nPsfParams; ++i) {
        _psf->getBasisKernel(i)->getImage(psfDerivative[i],  _center);
        Eigen::Map<VectoXd> basisImage = extractEigenView(psfDerivative[i]);
        basisImage *= _magnitude * coefficients[i];
    }

}
