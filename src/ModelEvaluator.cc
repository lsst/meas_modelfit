#include "lsst/meas/multifit/ModelEvaluator.h"
#include "lsst/meas/multifit/matrices.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/image/MaskedImage.h"

namespace multifit = lsst::meas::multifit;
namespace detection = lsst::afw::detection;

template <typename ImageT>
class CompressFunctor : 
    public detection::FootprintFunctor<lsst::afw::image::MaskedImage<ImageT> > {
public:
    typedef lsst::afw::image::MaskedImage<ImageT> MaskedImage;
    CompressFunctor(
        MaskedImage const & src,
        ndarray::Array<multifit::Pixel, 1, 1> const & imageDest,
        ndarray::Array<multifit::Pixel, 1, 1> const & varianceDest
    ) : detection::FootprintFunctor<MaskedImage>(src),
        _imageDest(imageDest),
        _varianceDest(varianceDest)        
    {}

    virtual void reset(detection::Footprint const & footprint) {
        if(_imageDest.getSize<0>() != footprint.getNpix() || 
            _varianceDest.getSize<0>() != footprint.getNpix()) {
            throw LSST_EXCEPT(lsst::pex::exceptions::LengthErrorException,
                "Destination vectors are not correct length");    
        }
        _imageIter = _imageDest.begin();
        _varianceIter = _varianceDest.begin(); 
    }
    virtual void operator()(
        typename MaskedImage::xy_locator loc,
        int x,
        int y
    ) {
       *_imageIter = loc.image();
       *_varianceIter = loc.variance();
       ++_imageIter;
       ++_varianceIter;
    }
private:
    ndarray::Array<multifit::Pixel, 1, 1> const & _imageDest;
    ndarray::Array<multifit::Pixel, 1, 1> const & _varianceDest;
    ndarray::Array<multifit::Pixel, 1, 1>::Iterator _imageIter, _varianceIter;
};

template class CompressFunctor<float>;
template class CompressFunctor<double>;


template <typename ImageT>
void multifit::ModelEvaluator::Traits<ImageT>::setExposureList(
    ModelEvaluator & evaluator,
    CalibratedExposureList const & exposureList
) {
    evaluator._projectionList.clear();
    evaluator._validProducts = 0;

    int nLinear = evaluator._model->getLinearParameterSize();
    int nNonlinear = evaluator._model->getNonlinearParameterSize();

    int pixSum;

    ModelProjection::Ptr projection;
    ExposureConstPtr exposure;
    FootprintConstPtr footprint;
    PsfConstPtr psf;
    WcsConstPtr wcs;    
  
    //exposures which contain fewer than _nMinPix pixels will be rejected
    //construct a list containing only those exposure which were not rejected
    std::list<ExposureConstPtr> goodExposureList;

    // loop to create projections
    for(typename CalibratedExposureList::const_iterator i(exposureList.begin()), 
        end(exposureList.end()); i != end; ++i
    ) {
        boost::tie(exposure, psf) = *i;
        wcs = exposure->getWcs();        
        footprint = evaluator._model->computeProjectionFootprint(
            psf, 
            wcs
        );
        MaskedImage const & maskedImage = exposure->getMaskedImage();
        footprint = fixFootprint(footprint, maskedImage.getMask());

        //ignore exposures with too few pixels        
        if (footprint->getNpix() >= evaluator._nMinPix) {
            ProjectionFrame frame(
                evaluator._model->makeProjection(psf, wcs, footprint)
            );
      
            evaluator._projectionList.push_back(frame);
            goodExposureList.push_back(exposure);

            pixSum += footprint->getNpix();
        }
    }

    //  allocate matrix buffers
    evaluator._imageVector = ndarray::allocate<Allocator>(ndarray::makeVector(pixSum));
    evaluator._varianceVector = ndarray::allocate<Allocator>(ndarray::makeVector(pixSum));

    evaluator._modelImage = ndarray::allocate<Allocator>(ndarray::makeVector(pixSum));
    evaluator._linearParameterDerivative = ndarray::allocate<Allocator>(
        ndarray::makeVector(nLinear, pixSum)
    );
    evaluator._nonlinearParameterDerivative = ndarray::allocate<Allocator>(
        ndarray::makeVector(nNonlinear, pixSum)
    );    
    
    int nPix;
    int pixelStart = 0, pixelEnd;

    typename std::list<ExposureConstPtr>::const_iterator exposureIter(
        goodExposureList.begin()
    );
    //loop to assign matrix buffers to each projection Frame
    for(ProjectionFrameList::iterator i(evaluator._projectionList.begin()), 
        end(evaluator._projectionList.end()); i != end; ++end
    ) {
        ProjectionFrame & frame(*i);
        nPix = frame.getFootprint()->getNpix();
        pixelEnd = pixelStart + nPix;

        // set image/variance buffers
        ndarray::shallow(frame._imageVector) = evaluator._imageVector[
            ndarray::view(pixelStart, pixelEnd)
        ];
        ndarray::shallow(frame._varianceVector) = evaluator._varianceVector[
            ndarray::view(pixelStart, pixelEnd)
        ];
        compressExposure(frame, *exposureIter); 
        
        //set modelImage buffer
        frame._projection->setModelImageBuffer(
            evaluator._modelImage[ndarray::view(pixelStart, pixelEnd)]
        );
        
        //set linear buffer
        frame._projection->setLinearParameterDerivativeBuffer(
            evaluator._linearParameterDerivative[ndarray::view()(pixelStart, pixelEnd)]
        );
        //set nonlinear buffer
        frame._projection->setNonlinearParameterDerivativeBuffer(
            evaluator._nonlinearParameterDerivative[ndarray::view()(pixelStart, pixelEnd)]
        );

        pixelStart = pixelEnd;
        ++exposureIter;
    }
}

template <typename ImageT>
void multifit::ModelEvaluator::Traits<ImageT>::compressExposure(
    ProjectionFrame & frame,
    ExposureConstPtr const & exposure 
) {
    CompressFunctor<ImageT> functor(
        exposure->getMaskedImage(), 
        frame._imageVector, 
        frame._varianceVector
    );
    functor.apply(*frame.getFootprint());    
}

template <typename ImageT>
multifit::FootprintConstPtr multifit::ModelEvaluator::Traits<ImageT>::fixFootprint(
    FootprintConstPtr const & footprint,
    typename MaskedImage::MaskPtr const & mask
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
        detection::Span const & span(**spanIter);
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
        typename Traits<ImageT>::MaskedImage::Mask::x_iterator maskIter = mask->x_at(x0, y);

        //loop over all span locations, slicing the span at maskedPixels
        for(int x = x0; x <= x1; ++x) {            
            if(*maskIter != 0) {
                //masked pixel found within span
                if (x > x0) {                    
                    //add beginning of spanto the fixedFootprint
                    //the fixed span contains all the unmasked pixels up to,
                    //but not including this masked pixel
                    fixedFootprint->addSpan(y, x0, x - 1);                
                }
                //set the next fixed Span to start after this pixel
                x0 = x + 1;
            }

            //move to next pixel
            ++maskIter;
        }
        //add last section of span
        if(x0 <= x1) {
            fixedFootprint->addSpan(y, x0, x1);
        }
    }
   
    fixedFootprint->setRegion(maskBBox);
    return boost::shared_ptr<Footprint const>(fixedFootprint);
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


