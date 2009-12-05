#include "lsst/meas/multifit/ModelEvaluator.h"
#include "lsst/meas/multifit/matrices.h"

namespace multifit = lsst::meas::multifit;
namespace detection = lsst::afw::detection;

class CompressFunctor : 
    detection::FootprintFunctor<multifit::MaskedImage> {
public:
    CompressFuntor(
        multifit::MaskedImage const & src,
        ndarray::Array<multifit::Pixel, 1, 1> const & imageDest,
        ndarray::Array<multifit::Pixel, 1, 1> const & varianceDest
    ) : lsst::afw::detection::FootprintFunctor(src),
        _imageDest(imageDest),
        _varianceDest(varianceDest),        
    {}

    virtual void reset(detection::Footprint const & footprint) {
        if(imageDest.getSize() != footprint.getNpix() || 
            varianceDest.getSize() != footprint.getNPix()) {
            throw LSST_EXCEPT(lsst::pex::exceptions::LengthErrorException,
                "Destination vectors are not correct length");    
        }
        _imageIter = _imageDest.begin();
        _variaceIter = _varianceDest.begin(); 
    }
    virtual void operator()(
        typename multifit::ModelEvaluator::MaskedImage::xy_locator loc,
        int x,
        int y
    ) {
       *_imageIter = loc.image();
       *_varianceIter = loc.var();
       ++_imageIter;
       ++_varianceIter;
    }
private:
    ndarray::Array<multifit::Pixel, 1, 1> const & _imageDest;
    ndarray::Array<multifit::Pixel, 1, 1> const & _varianceDest;
    ndarray::Array<multifit::Pixel, 1, 1>::iterator _imageIter, _varianceIter;


}

void multifit::ModelEvaluator::Frame::compressExposure() {
    if(!exposure || !imageVector || !variaceVector) {
        return;
    }
    CompressFunctor functor(
        *exposure.getMaskedImage(), imageVector, varianceVector
    );
    functor.apply(*footprint);    
}

void multifit::ModelEvaluator::setExposureList(
    CalibratedExposureList const & exposureList,
    int nMinPix
) {
    int nExposures = exposureList.size();
    _projectionList.clear();


    int nLinear = _model->getLinearParameterSize();
    int nNonLinear = _model->getNonlinearParameterSize();

    int pixSum;

    ModelProjectionPtr projection;
    ExposurePtr exposurePtr;
    MaskedImagePtr maskedImagePtr;
    FootprintPtr footprintPtr;
    KernelPtr kernelPtr;
    WcsPtr wcsPtr;    
    double photFactor;    
    
    CalibratedExposureIterator exposureIter(exposuresList.begin());
    CalibratedExposureIterator const & exposureEnd(exposureList.end());

    // loop to create projections
    for( ; exposureIter != exposureEnd; ++exposureIter) {
        boost::tie(exposurePtr, kernelPtr, photFactor) = *exposureIter;
        wcsPtr = exposurePtr->getWcs();        
        footprintPtr = _model->computeProjectionFootprint(
            *kernelPtr, 
            wcsPtr,
            photoFactor
        );
        maskedImagePtr = exposurePtr->getMaskedImage();
        footprintPtr = fixFootprint(footprintPtr, maskedImagePtr);

        //ignore exposures with too few pixels        
        if (footprintPtr->getNPix() < nMinPix) {
            continue
        }

        ModelProjectionPtr projectionPtr = _model->makeProjection(
            *kernelPtr, wcsPtr, footprintPtr, photFactor, _activeProducts
        );

        Frame frame(projectionPtr, exposurePtr, kernelPtr, 
            wcsPtr, footprintPtr);
      
        _projectionList.push_back(frame);
        pixSum += footprintPtr->getNpix();
    }

    //  allocate matrix buffers
    _imageVector = ndarray::allocate<Allocator>(pixSum);
    _varianceVector = ndarray::allocate<Allocator>(pixSum);

    if(_activeProducts & ModelProjection::MODEL_IMAGE) {
        _modelImage = ndarray::allocate<Allocator>(pixSum);
    }
    if(_activeProducts & ModelProjection::LINEAR_PARAMETER_DERIVATIVE) {
        _linearParameterDerivative = ndarray::allocate<Allocator>(
            ndarray::MakeVector(nLiear, pixSum)
        );
    }
    if(_activeProducts & ModelProjection::NONLINEAR_PARAMETER_DERIVATIVE) {
        _nonlinearParameterDerivative = ndarray::allocate<Allocator>(
            ndarray::MakeVector(nNonlinear, pixSum)
        );    
    }
    
    int nPix;
    int pixelStart = 0, pixelEnd;

    FrameIterator frameIter(_projectionList.begin());
    FrameIterator const & FrameEnd(_projectionList.end());
        
    //loop to assign matrix buffers to each projection Frame
    for( ; frameIter != frameEnd; ++frameIter) {
        Frame & frame(*frameIter);
        nPip = frame.footprint->getNpix();
        pixelEnd = pixelStart + nPix;

        frame.imageVector = _imageVector[ndarray::view(pixelStart, pixelEnd)];
        frame.varianceVector = 
            _varianceVector[ndarray::view(pixelStart, pixelEnd)];
        frame.compressExposure();

        if(_activeProducts & ModelProjection::MODEL_IMAGE) {            
            frame.modelImage = _modelImage[ndarray::view(pixelStart, pixelEnd)];
            frame.projection->setModelImageBuffer(frame.modelImage);            
        }
        if(_activeProducts & ModelProjection::LINEAR_PARAMETER_DERIVATIVE) {
            frame.linearParameterDerivative = _linearParameterDerivative[
                ndarray::view()(pixelStart, pixelEnd)
            ];
            frame.projection->setLinearParameterDerivativeBuffer(
                frame.linearParameterDerivative
            );
        }
        if(_activeProducts & ModelProjection::NONLINEAR_PARAMETER_DERIVATIVE) {
            frame.nonlinearParameterDerivative = _nonlinearParameterDerivative[
                ndarray::view()(pixelStart, pixelEnd)
            ];
            frame.projection->setNonlinearParameterDerivativeBuffer(
                frame.nonlinearParameterDerivative
            );
        }
        if(_activeProducts & ModelProjection::PSF_PARAMETER_DERIVATIVE) {
            int nPsf = frame.projection->getPsfParameterSize();
            frame.linearParameterDerivative = ndarray::allocate<Allocator>(
                ndarray::makeVector(nPsf, nPix)
            );
            frame.projection->setPsfParameterDerivativeBuffer(
                frame.psfParameterDerivative
            );
        }
        if(_activeProducts & ModelProjection::WCS_PARAMETER_DERIVATIVE) {
            int nWcs = frame.projection->getWcsParameterSize();
            frame.wcsParameterDerivative = ndarray::allocate<Allocator> (
                ndarray::makeVector(nWcs, nPix)
            );
            frame.projection->setWcsParameterDerivativeBuffer(
                frame.wcsParameterDerivative
            );
        }

        pixelStart = pixelEnd;
    }

}

multifit::ModelEvaluator::FootprintPtr multifit::ModelEvaluator::fixFootprint(
    FootprintPtr const & footprint,
    MaskPtr const & mask
) {
    lsst::afw::image::BBox const maskBBox(
        mask->getXY0(), mask->getWidth(), mask->getHeight()
    );
    int maskX0 = maskBBox.getX0();
    int maskY0 = maskBBox.getY0();
    int maskX1 = maskBBox.getX1();
    int maskY1 = maskBBox.getY1();
     
    Footprint * fixedFootprint = new Footprint();    
    Footprint::SpanList const & oldSpanList(footprint->getSpans());

    Footprint::SpanList::const_iterator spanIter(oldSpanList.begin());
    Footprint::SpanList::const_iterator const & spanEnd(oldSpanList.end());
    

    int x0, x1, y;
    for( ; spanIter != spanEnd; ++spanIter) {
        Footprint::Span const & span(**spanIter);
        y = span.getY();
        x0 = span.getX0();
        x1 = span.getX1();
        if(y < maskY0 || y > maskY1 || x1 < maskX0 || x0 > maskX1) {
            //span is entirely outside the image mask. 
            //cannot be used
            continue;
        }
        if(x0 < maskX0) x0 = maskX0;
        if(x1 > maskX1) x1 = maskX1;
        Mask::x_iterator maskIter = 
            mask.x_at(spanIter->getX0(), spanIter->getY());

        //loop over all span locations, slicing the span at maskedPixels
        for(int x = x0; x <= x1; ++x) {            
            if(*maskIter != 0) {
                //masked pixel found within span
                if (x > x0) {                    
                    //add beginning of spanto the fixedFootprint
                    //the fixed span contains all the unmasked pixels up to,
                    //but not including this masked pixel
                    fixedFootprint.addSpan(y, x0, x - 1);                
                }
                //set the next fixed Span to start after this pixel
                x0 = x + 1;
            }

            //move to next pixel
            ++maskIter;
        }
        //add last section of span
        if(x0 <= x1) {
            fixedFootprint.addSpan(y, x0, x1);
        }
    }
   
    fixedFootprint->setRegion(maskBBox);
    return boost::shared_ptr<Footprint const>(fixedFootprint);
}

ndarray::Array<multifit::Pixel const, 1, 1> multifit::ModelEvaluator::computeModelImage() {
    if(!(_activeProducts & ModelProjection::MODEL_IMAGE)) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Product MODEL_IMAGE is not enabled."
        );
    }
    if(!(_validProducts & ModelProjection::MODEL_IMAGE)) {
        FrameIterator i(_projectionList.begin());
        FrameIterator const & end(_projectionList.end());
        for( ; i  != end; ++i) {
            i->projection->computeModelImage();
        }
    }    
    return _modelImage;
}

ndarray::Array<multifit::Pixel const, 2, 2> multifit::ModelEvaluator::computeLinearParameterDerivative() {
    if(!(_activeProducts & ModelProjection::LINEAR_PARAMETER_DERIVATIVE)) {     
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Product MODEL_IMAGE is not enabled."
        );
    }
    if (!(_validProducts & LINEAR_PARAMETER_DERIVATIVE)) {
        FrameIterator i(_projectionList.begin());
        FrameIterator const & end(_projectionList.end());
        for( ; i  != end; ++i) {
            i->projection->computeLinearParameterDerivative();
        }
    }    
    return _linearParameterDerivative;
}

ndarray::Array<multifit::Pixel const, 2, 2> multifit::ModelEvaluator::computeNonlinearParameterDerivative() {
    if (!(_activeProducts & NONLINEAR_PARAMETER_DERIVATIVE)) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Product NONLINEAR_PARAMETER_DERIVATIVE is not enabled."
        );
    }
    if(!(_validProducts & ModelProjection::NONLINEAR_PARAMETER_DERIVATIVE)) {
        FrameIterator i(_projectionList.begin());
        FrameIterator const & end(_projectionList.end());
        for( ; i  != end; ++i) {
            i->projection->computeNonlinearParameterDerivative();
        }
    }    
    return _nonlinearParameterDerivative;
}

ndarray::Array<multifit::Pixel const, 2, 2> multifit::ModelEvaluator::computeWcsParameterDerivative() {
    if (!(_activeProducts & WCS_PARAMETER_DERIVATIVE)) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Product WCS_PARAMETER_DERIVATIVE is not enabled."
        );
    }
    if(!(_validProducts & ModelProjection::WCS_PARAMETER_DERIVATIVE)) {
        FrameIterator i(_projectionList.begin());
        FrameIterator const & end(_projectionList.end());
        for( ; i  != end; ++i) {
            i->projection->computeWcsParameterDerivative();
        }
    }    
    return _wcsParameterDerivative;
}

ndarray::Array<multifit::Pixel const, 2, 2> multifit::ModelEvaluator::computePsfParameterDerivative() {
    if(!(_activeProducts & ModelProjection::PSF_PARAMETER_DERIVATIVE)) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Product PSF_PARAMETER_DERIVATIVE is not enabled."
        );    
    }
    if(!(_validProducts & ModelProjection::PSF_PARAMETER_DERIVATIVE)) {
        FrameIterator i(_projectionList.begin());
        FrameIterator const & end(_projectionList.end());
        for( ; i  != end; ++i) {
            i->projection->computePsfParameterDerivative();
        }
    }    
    return _linearParameterDerivative;
}

#if 0 
//old code. kept for reference

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
#endif
