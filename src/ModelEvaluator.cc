#include "lsst/meas/multifit/ModelEvaluator.h"
#include "lsst/meas/multifit/matrices.h"
#include "lsst/meas/multifit/footprintUtils.h"
#include "lsst/afw/image/MaskedImage.h"

#include <iostream>
namespace multifit = lsst::meas::multifit;

template<typename ImagePixel, typename MaskPixel, typename VariancePixel>
void multifit::ModelEvaluator::setExposureList(
    std::list< boost::shared_ptr< CharacterizedExposure<
            ImagePixel, MaskPixel, VariancePixel
        > > > const & exposureList
) { 
    typedef std::list< 
        boost::shared_ptr< 
            CharacterizedExposure<ImagePixel,MaskPixel,VariancePixel> 
        > 
    > ExposureList;
    typedef typename ExposureList::const_iterator ExposureIterator;
    
    _projectionList.clear();
    _validProducts = 0;
    
    int nLinear = getLinearParameterSize();
    int nNonlinear = getNonlinearParameterSize();

    int pixSum = 0;

    typename CharacterizedExposure<ImagePixel,MaskPixel,VariancePixel>::Ptr exposure;
    ModelProjection::Ptr projection;
    FootprintConstPtr footprint;
    PsfConstPtr psf;
    WcsConstPtr wcs;    
  
    //exposures which contain fewer than _nMinPix pixels will be rejected
    //construct a list containing only those exposure which were not rejected
    ExposureList goodExposureList;

    // loop to create projections
    for(ExposureIterator i(exposureList.begin()), end(exposureList.end());
        i != end; ++i
    ) {
        exposure = *i;
        psf = exposure->getPSF();
        wcs = exposure->getWcs();        
        footprint = _model->computeProjectionFootprint(psf, wcs);

        footprint = clipAndMaskFootprint<MaskPixel>(
            *footprint, exposure->getMaskedImage().getMask()
        );
        //ignore exposures with too few contributing pixels        
        if (footprint->getNpix() > _nMinPix) {
            ProjectionFrame frame(
                _model->makeProjection(psf, wcs, footprint)
            );
      
            _projectionList.push_back(frame);
            goodExposureList.push_back(exposure);

            pixSum += footprint->getNpix();
        }
    }

    //  allocate matrix buffers
    ndarray::shallow(_imageVector) = ndarray::allocate<Allocator>(ndarray::makeVector(pixSum));
    ndarray::shallow(_varianceVector) = ndarray::allocate<Allocator>(ndarray::makeVector(pixSum));

    ndarray::shallow(_modelImage) = ndarray::allocate<Allocator>(ndarray::makeVector(pixSum));
    ndarray::shallow(_linearParameterDerivative) = ndarray::allocate<Allocator>(
        ndarray::makeVector(nLinear, pixSum)
    );
    ndarray::shallow(_nonlinearParameterDerivative) = ndarray::allocate<Allocator>(
        ndarray::makeVector(nNonlinear, pixSum)
    );    
    
    int nPix;
    int pixelStart = 0, pixelEnd;

    ExposureIterator exposureIter(goodExposureList.begin());

    //loop to assign matrix buffers to each projection Frame
    for(ProjectionFrameList::iterator i(_projectionList.begin()), end(_projectionList.end()); 
        i != end; ++i
    ) {
        ProjectionFrame & frame(*i);
        nPix = frame.getFootprint()->getNpix();
        pixelEnd = pixelStart + nPix;

        // set image/variance buffers
        ndarray::shallow(frame._imageVector) = _imageVector[
            ndarray::view(pixelStart, pixelEnd)
        ];
        ndarray::shallow(frame._varianceVector) = _varianceVector[
            ndarray::view(pixelStart, pixelEnd)
        ];

        // compress the exposure using the footprint
        compressImage(
            *frame.getFootprint(), 
            (*exposureIter)->getMaskedImage(), 
            frame._imageVector, 
            frame._varianceVector
        );

        //set modelImage buffer
        frame._projection->setModelImageBuffer(
            _modelImage[ndarray::view(pixelStart, pixelEnd)]
        );
        
        //set linear buffer
        frame._projection->setLinearParameterDerivativeBuffer(
            _linearParameterDerivative[ndarray::view()(pixelStart, pixelEnd)]
        );
        //set nonlinear buffer
        frame._projection->setNonlinearParameterDerivativeBuffer(
            _nonlinearParameterDerivative[ndarray::view()(pixelStart, pixelEnd)]
        );

        
        pixelStart = pixelEnd;
        ++exposureIter;
    }   
}

ndarray::Array<multifit::Pixel const, 1, 1> multifit::ModelEvaluator::computeModelImage() {
    if(!(_validProducts & MODEL_IMAGE)) {
        ProjectionFrameIterator i(_projectionList.begin());
        ProjectionFrameIterator const & end(_projectionList.end());
        for( ; i  != end; ++i) {
            i->computeModelImage();
        }
    }    
    return _modelImage;
}

ndarray::Array<multifit::Pixel const, 2, 2> multifit::ModelEvaluator::computeLinearParameterDerivative() {
    if (!(_validProducts & LINEAR_PARAMETER_DERIVATIVE)) {
        ProjectionFrameIterator i(_projectionList.begin());
        ProjectionFrameIterator const & end(_projectionList.end());
        for( ; i  != end; ++i) {
            i->computeLinearParameterDerivative();
        }
        _validProducts |= LINEAR_PARAMETER_DERIVATIVE;
    }    
    return _linearParameterDerivative;
}

ndarray::Array<multifit::Pixel const, 2, 2> multifit::ModelEvaluator::computeNonlinearParameterDerivative() {
    if(!(_validProducts & NONLINEAR_PARAMETER_DERIVATIVE)) {
        ProjectionFrameIterator i(_projectionList.begin());
        ProjectionFrameIterator const & end(_projectionList.end());
        for( ; i  != end; ++i) {
            i->computeNonlinearParameterDerivative();
        }
        _validProducts |= NONLINEAR_PARAMETER_DERIVATIVE;
    }    
    return _nonlinearParameterDerivative;
}


template void multifit::ModelEvaluator::setExposureList<float, 
    lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>
(
    std::list<boost::shared_ptr<CharacterizedExposure<float> > > const &
);
template void multifit::ModelEvaluator::setExposureList<double, 
    lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>
(
    std::list<boost::shared_ptr<CharacterizedExposure<double> > > const &
);
