#include "lsst/meas/multifit/ModelEvaluator.h"

namespace multifit = lsst::meas::multifit;


template <typename ExposureContainer>
multifit::ModelEvaluator::ModelEvaluator(
        ObjectModel::Map const & model_set, 
        ExposureContainer exposures,
        int marginalizationFlags,
        Probability const * prior) {
    //TODO
}

template <typename ExposureContainer>
multifit::ModelEvaluator::ModelEvaluator(
        ObjectModel::Map const & model_set,
        ndarray::ArrayRef<double,1,1> const & linear,
        ndarray::ArrayRef<double,1,1> const & nonlinear,
        ExposureContainer exposures,
        int marginalizationFlags,
        Probability const * prior) :
    _models(model_set),
    _linear(linear.core().vector()), 
    _nonlinear(nonlinear.core().vector()),
    _marginalizationFlags(marginalizationFlags) 
{
    setExposures(exposures);

    //build background subtracted data
    buildBkgSubtracted()
    
    //build list of SourceModels.
    //this is a kind of workspace for matrix calculations
    _workspace.reserve(_exposures.size()*_models.size());
    
    ExposureIndexIterator exposure_iter(_exposures.begin());
    ExposureIndexIterator const exposure_end(_exposures.end());
    ObjectModelIterator const object_end(_models.end());
    ObjectModelIterator object_iter(_models.begin());
    exposure_iter = _exposures.begin();
    for(; exposure_iter != exposure_end; ++exposure_iter) {
        object_iter = _models.begin();
        for(;object_iter != object_end; ++object_iter) {
            SourceModel * source_model = (*object_iter)->createSource(
                    exposure_iter->first, 
                    _marginalizationFlags);

            _workspace.push_back(source_model);
        }
    }
    

    if(prior != NULL) {
        Eigen::VectorXd params(_nonlinear.size() + _linear.size());
        params << _nonlinear, _linear;

        _posterior = prior->at(params);
    }
 }

template <typename ExposureContainer>
void multifit::ModelEvaluator::setExposures(
        ExposureContainer const & exposures) {
    _exposures.clear();
    _residualsSection.clear();
    _linearMatrixSection.clear();
    _nonlinearMatrixSection.clear();
    _psfMatrixSection.clear();

    int nExposures = exposures.size();
    _exposures.reserve(nExposures);
    _residualsSection.reserve(nExposures);
    _linearMatrixSection.reserve(nExposures);
    _nonlinearMatrixSection.reserve(nExposures);
    _psfMatrixSection.reserve(nExposures);

    _numTotalPixels=0;
    int nPix;
    ExposureContainer::const_iterator iter = exposures.begin();
    ExposureContainer::const_iterator const end = exposures.end();
    for(; iter != end; ++iter) {
        CalibratedExposure::Ptr exposure = *iter;
        _exposures.push_back(exposure);
        nPix = exposure->getWidth()*exposure->getHeight();
        _bkgSubtracted.segment(_numTotalPixels, nPix) = 
                iter->getSigmaWeightedBackgroundSubtracted();
        _numTotalPixels += nPix;
    }

    _linearMatrix.resize(_linear.size(), _numTotalPixels);
    _nonlinearMatrix.resize(_nonlinear.size(), _numTotalPixels);
    _residuals.resize(_numTotalPixels);
    _bkgSubtracted.resize(_numTotalPixels);

    int linearMatrixOffset = 0;
    int nonLinearMatrixOffset = 0;
    int psfMatrixOffset = 0;
    int residualsOffset = 0;
    int height, width;

    for(iter = exposures.begin(); iter != end; ++iter) {
        height = iter->getHeight();
        width = iter->getWidth();
        nPix = height*width;
        
        _residualsSection.push_back(ImageCore(
                _residuals.data() + residualsOffset,
                ndarray::make_index(height, width),
                ndarray::make_index(width, 1))
        );
        residualsOffset += nPix;

        _linearMatrixSection.push_back(DerivativeCore(
                _linearMatrix.data() + linearMatrixOffset,
                ndarray::make_index(_linear.size(), height, width),
                ndarray::make_index(nPix, width, 1))
        );
        linearMatrixOffset +=_linear.size()*nPix;

        _nonlinearMatrixSection.push_back(DerivativeCore(
                _nonlinearMatrix.data() + nonlinearMatrixOffset,
                ndarray::make_index(_nonlinear.size(), height, width),
                ndarray::make_index(nPix, width, 1))
        );
        nonlinearMatrixOffset +=_linear.size()*nPix;
    }
}

void multifit::ModelEvaluator::buildBkgSubtracted() {
    //TODO: implement me
    //may require doing background subtraction
    //then copying into eigen vector
}

multifit::ProbabilityExpansion const & multifit::ModelEvaluator::evaluate() {
    int const linearSize = _linear.size();
    int const nonlinearSize = _nonlinear.size();
    int const nExposures = _exposures.size();

    _residuals = _bkgSubtracted;

    ParameterVector linearRef(_linear.data(), _linear.size(), 1);
    ParameterVector nonlinearRef(_nonlinear.data(), _nonlinear.size(), 1);

    //loop to compute linear matrix    

    for(int exposureId=0; exposureId < nExposures; ++exposureId) {
        TransformedModel::Ptr model = _modelStack.at(exposureId);
        model->evalLinearDerivative(
                nonlinearRef, 
                getLinearMatrixSection(exposureId),
                getResidualsSection(exposureId)
        );
    }

    //solve for linear params
    Eigen::Block<Eigen::MatrixXd> ddpLinearBlock(_posterior._ddp, 
            nonlinearSize, 
            nonlinearSize, 
            linearSize, 
            linearSize);

    ddpLinearBlock.part<Eigen::SelfAdjoint>() += 
            (_linearMatrix.adjoint() * _linearMatrix).lazy();

    //TODO: implement solver
    //solve(ddpLinearBlock, 
    //        _linearMatrix.transpose()*_residuals - _posterior._dp.end(linearSize), 
    //        _linear);
                    
    //compute residuals
    _residuals -= _linearMatrix * _linear;
    
    //loop to compute nonlinear matrix
    source_iter = _workspace.begin();
    exposure_iter = _exposures.begin();
    for(int exposureId = 0; exposureId < nExposures; ++exposureId) {
        TransformedModel::Ptr model = _modelStack.at(exposureId);
        model->evalNonlinearDerivative(
                linearRef,
                nonlinearRef,
                getNonlinearMatrixSection(exposureId)
        );                        
    }

    _posterior._p += 0.5*_residuals.dot(_residuals);
    _posterior._dp.end(linearSize).setZero();
        _posterior._dp.start(nonlinearSize) -= 
            _nonlinearMatrix.transpose() * _residuals;
    _posterior._ddp.block(0, 0, nonlinearSize, nonlinearSize).part<Eigen::SelfAdjoint>() += 
            (_nonlinearMatrix.adjoint() * _nonlinearMatrix).lazy();
    _posterior._ddp.block(nonlinearSize, 0, linearSize, nonlinearSize) +=
            _linearMatrix.adjoint() * _nonlinearMatrix;

    //loop to marginalize over calibration parameters...lots of work to do
    source_iter = _workspace.begin();
    exposure_iter = _exposures.begin();
    for(int exposureId ; exposureId < nExposures; ++exposureId) {
        int width = exposure_iter->first->getWidth();
        int height = exposure_iter->first->getHeight();
        int pixels = width*height;
        int bkg_size = getNumBkgParam(exposure_iter->first);
        int psf_size = getNumPsfParam(exposure_iter->first);
        int calibration_size = bkg_size + psf_size;
        if (calibration_size > 0) {
            _calibration_matrix.resize(width*height,calibration_size);
            _calibration_matrix.setZero();
        }
        object_iter = ObjectModelIterator::begin(_models);
        for(; object_iter != object_end; ++object_iter, ++source_iter) {
            if (psf_size > 0) {  
                (*object_iter)->computePsfMatrix(
                        *source_iter, 
                        getPsfMatrix(exposure_iter)
                );
            }
        }
        if (calibration_size > 0) {
            if (bkg_size > 0) computeBkgMatrix(exposure_iter->first);

            Eigen::VectorXd dtau = _calibration_matrix.transpose() * _residuals;
            Eigen::MatrixXd dtau_dtau(calibration_size,calibration_size);
            dtau_dtau.part<Eigen::SelfAdjoint>() 
                = (_calibration_matrix.adjoint() * _calibration_matrix).lazy();
            Eigen::MatrixXd dtau_dphi = _calibration_matrix.adjoint() 
                    * _nonlinearMatrix.block(
                            exposure_iter->second,0,pixels,nonlinearSize);
            Eigen::MatrixXd dtau_dlambda = _calibration_matrix.adjoint() 
                    * _linearMatrix.block(
                            exposure_iter->second,0,pixels,linearSize);
            addCalibrationFisherMatrix(exposure_iter->first,dtau_dtau);
            Eigen::LU<Eigen::MatrixXd> lu(dtau_dtau.lu());
            _posterior._p += 0.5*std::log(
                    lu.determinant() * std::pow(2*M_PI,-calibration_size)
            );
            
            Eigen::MatrixXd dtau_dtau_inv;

            // TODO: can optimize stuff below by using Cholesky instead of lU, 
            // storing some temporaries.
            lu.computeInverse(&dtau_dtau_inv);

            double temp = (dtau.adjoint() * dtau_dtau_inv).dot(dtau);
            _posterior._p -= 0.5 * temp;
            _posterior._dp.start(nonlinearSize) += dtau_dphi.adjoint() * dtau_dtau_inv * dtau;
            _posterior._dp.end(linearSize) += dtau_dlambda.adjoint() * dtau_dtau_inv * dtau;
            Eigen::Block<Eigen::MatrixXd> ddp_nonlinear = 
                    _posterior._ddp.block(0, 0, nonlinearSize, nonlinearSize);
            ddp_nonlinear.part<Eigen::SelfAdjoint>() += 
                    (dtau_dphi.adjoint()*dtau_dtau_inv.part<Eigen::SelfAdjoint>()*dtau_dphi).lazy();

            Eigen::Block<Eigen::MatrixXd> ddp_linear = _posterior._ddp.block(
                    nonlinearSize, nonlinearSize, linearSize,linearSize);
            ddp_linear.part<Eigen::SelfAdjoint>() +=
                (dtau_dlambda.adjoint()*dtau_dtau_inv.part<Eigen::SelfAdjoint>()*dtau_dlambda).lazy();

            _posterior._ddp.block(nonlinearSize,0,linearSize,nonlinearSize) +=
                dtau_dlambda.adjoint() * dtau_dtau_inv * dtau_dphi;
            
        }
    }

    // make 
    _posterior._ddp.block(0, nonlinearSize, nonlinearSize, linearSize) =
        _posterior._ddp.block(nonlinearSize, 0, linearSize, nonlinearSize).adjoint();

    
    return _posterior;
}

/**
 */
Eigen::VectorXd multifit::ModelEvaluator::step() {
    Eigen::VectorXd params(_linear.size() + _nonlinear.size());

    //solve(_posterior.ddp, -_posterior.dp, params);
    _nonlinear += params.start(_nonlinear.size());
    _linear += params.end(_linear.size());
    return params;
        
}
