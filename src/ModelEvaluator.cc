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
    generateExposureIndex(exposures.begin(), exposures.end());
    _linear_matrix.resize(_linear.size(), _numTotalPixels);
    _nonlinear_matrix.resize(_nonlinear.size(), _numTotalPixels);
    _residuals.resize(_numTotalPixels);
    _bkg_subtracted.resize(_numTotalPixels);
    
    //build background subtracted data
    buildBkgSubtracted();
    
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

template <typename ExposureIterator>
void multifit::ModelEvaluator::setExposures(
        ExposureIterator const & begin,
        ExposureIterator const & end) {
    _exposures.empty();

    int offset = 0;
    int size;
    for(ExposureIterator iter = begin; iter != end; iter++) {
        _exposures.insert(std::make_pair(*iter, offset));        
        size = iter->getWidth()*iter->getHeight();
        _bkg_subtracted.segment(offset, size) = 
                iter->getSigmaWeightedBackgroundSubtracted();
        offset += size;
    }
    _numTotalPixels = offset;
}

void multifit::ModelEvaluator::buildBkgSubtracted() {
    //TODO: implement me
    //may require doing background subtraction
    //then copying into eigen vector
}

multifit::ProbabilityExpansion const & multifit::ModelEvaluator::evaluate() {
    const int linear_size = _linear.size();
    const int nonlinear_size = _nonlinear.size();

    _residuals = _bkg_subtracted;

    ExposureIndexIterator exposure_iter(_exposures.begin());
    ExposureIndexIterator const exposure_end(_exposures.end());
    ObjectModelIterator object_iter;
    ObjectModelIterator const object_end(_models.end());
    SourceModelIterator source_iter(_workspace.begin());    
    //loop to compute linear matrix        
    for( ; exposure_iter != exposure_end; ++exposure_iter) {
        object_iter = ObjectModelIterator::begin(_models);
        for( ; object_iter != object_end; ++object_iter, ++source_iter) { 
            ObjectModel::Ptr model = *object_iter;
            model->computeLinearMatrix(**source_iter, 
                    getNonlinearParameterSection(object_iter),
                    getLinearMatrixSection(object_iter, exposure_iter),
                    getResidualsSection(exposure_iter));
        }        
    }

    //solve for linear params
    Eigen::Block<Eigen::MatrixXd> ddp_linear_block(_posterior._ddp, 
            nonlinear_size, 
            nonlinear_size, 
            linear_size, 
            linear_size);

    ddp_linear_block.part<Eigen::SelfAdjoint>() += 
            (_linear_matrix.adjoint() * _linear_matrix).lazy();

    //TODO: implement solver
    //solve(ddp_linear_block, 
    //        _linear_matrix.transpose()*_residuals - _posterior._dp.end(linear_size), 
    //        _linear);
                    
    //compute residuals
    _residuals -= _linear_matrix * _linear;
    
    //loop to compute nonlinear matrix
    source_iter = _workspace.begin();
    exposure_iter = _exposures.begin();
    for(; exposure_iter != exposure_end; ++exposure_iter) {
        object_iter = ObjectModelIterator::begin(_models);
        for(; object_iter != object_end; ++object_iter) {
            ObjectModel::Ptr model = *object_iter;
            model->computeNonlinearMatrix(**source_iter,
                    getLinearParameterSection(object_iter),
                    getNonlinearMatrixSection(object_iter, exposure_iter));
            ++source_iter;
        }
    }

    _posterior._p += 0.5*_residuals.dot(_residuals);
    _posterior._dp.end(linear_size).setZero();
    _posterior._dp.start(nonlinear_size) -= 
            _nonlinear_matrix.transpose() * _residuals;
    _posterior._ddp.block(0, 0, nonlinear_size, nonlinear_size).part<Eigen::SelfAdjoint>() += 
            (_nonlinear_matrix.adjoint() * _nonlinear_matrix).lazy();
    _posterior._ddp.block(nonlinear_size, 0, linear_size, nonlinear_size) +=
            _linear_matrix.adjoint() * _nonlinear_matrix;

    //loop to marginalize over calibration parameters...lots of work to do
    source_iter = _workspace.begin();
    exposure_iter = _exposures.begin();
    for( ; exposure_iter != exposure_end; ++exposure_iter) {
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
                    * _nonlinear_matrix.block(
                            exposure_iter->second,0,pixels,nonlinear_size);
            Eigen::MatrixXd dtau_dlambda = _calibration_matrix.adjoint() 
                    * _linear_matrix.block(
                            exposure_iter->second,0,pixels,linear_size);
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
            _posterior._dp.start(nonlinear_size) += dtau_dphi.adjoint() * dtau_dtau_inv * dtau;
            _posterior._dp.end(linear_size) += dtau_dlambda.adjoint() * dtau_dtau_inv * dtau;
            Eigen::Block<Eigen::MatrixXd> ddp_nonlinear = 
                    _posterior._ddp.block(0, 0, nonlinear_size, nonlinear_size);
            ddp_nonlinear.part<Eigen::SelfAdjoint>() += 
                    (dtau_dphi.adjoint()*dtau_dtau_inv.part<Eigen::SelfAdjoint>()*dtau_dphi).lazy();

            Eigen::Block<Eigen::MatrixXd> ddp_linear = _posterior._ddp.block(
                    nonlinear_size, nonlinear_size, linear_size,linear_size);
            ddp_linear.part<Eigen::SelfAdjoint>() +=
                (dtau_dlambda.adjoint()*dtau_dtau_inv.part<Eigen::SelfAdjoint>()*dtau_dlambda).lazy();

            _posterior._ddp.block(nonlinear_size,0,linear_size,nonlinear_size) +=
                dtau_dlambda.adjoint() * dtau_dtau_inv * dtau_dphi;
            
        }
    }

    // make 
    _posterior._ddp.block(0, nonlinear_size, nonlinear_size, linear_size) =
        _posterior._ddp.block(nonlinear_size, 0, linear_size, nonlinear_size).adjoint();

    
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
