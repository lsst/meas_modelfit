#ifndef LSST_MEAS_MULTIFIT_FOOTPRINT_UTILS_H
#define LSST_MEAS_MULTIFIT_FOOTPRINT_UTILS_H

#include "boost/shared_ptr.hpp"
#include "ndarray.hpp"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/meas/multifit/core.h"

namespace lsst {
namespace meas {
namespace multifit {

lsst::afw::detection::Footprint::Ptr makeFootprint(
    lsst::afw::geom::ellipses::Ellipse const & ellipse
);

template <typename MaskPixel>
lsst::afw::detection::Footprint::Ptr clipAndMaskFootprint(
    FootprintConstPtr const & footprint,
    typename lsst::afw::image::Mask<MaskPixel>::Ptr const & mask
) {
    lsst::afw::image::BBox const maskBBox(
        mask->getXY0(), mask->getWidth(), mask->getHeight()
    );
    int maskX0 = maskBBox.getX0();
    int maskY0 = maskBBox.getY0();
    int maskX1 = maskBBox.getX1();
    int maskY1 = maskBBox.getY1();
     
    lsst::afw::detection::Footprint * fixedFootprint = new Footprint();    
    lsst::afw::detection::Footprint::SpanList const & oldSpanList(
        footprint->getSpans()
    );

    lsst::afw::detection::Footprint::SpanList::const_iterator spanIter(
        oldSpanList.begin()
    );
    lsst::afw::detection::Footprint::SpanList::const_iterator const & spanEnd(
        oldSpanList.end()
    );
    

    int x0, x1, y;
    for( ; spanIter != spanEnd; ++spanIter) {
        lsst::afw::detection::Span const & span(**spanIter);
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
        typename lsst::afw::image::Mask<MaskPixel>::x_iterator maskIter = mask->x_at(x0, y);

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
    return boost::shared_ptr<lsst::afw::detection::Footprint>(fixedFootprint);
}


template <typename ImagePixel, typename MaskPixel=lsst::afw::image::MaskPixel,
    typename VariancePixel=lsst::afw::image::VariancePixel>
class CompressFunctor : 
    public lsst::afw::detection::FootprintFunctor<lsst::afw::image::MaskedImage<ImagePixel, MaskPixel, VariancePixel> > {
public:
    typedef lsst::afw::image::MaskedImage<ImagePixel, MaskPixel, VariancePixel> MaskedImage;
    typedef lsst::afw::detection::FootprintFunctor<MaskedImage> FootprintFunctor;

    CompressFunctor(
        MaskedImage const & src,
        ndarray::Array<Pixel, 1, 1> const & imageDest,
        ndarray::Array<Pixel, 1, 1> const & varianceDest
    ) : FootprintFunctor(src),
        _imageDest(imageDest),
        _varianceDest(varianceDest)        
    {}

    virtual void reset(lsst::afw::detection::Footprint const & footprint) {
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

    ndarray::Array<Pixel, 1, 1> const & _imageDest;
    ndarray::Array<Pixel, 1, 1> const & _varianceDest;
    ndarray::Array<Pixel, 1, 1>::Iterator _imageIter, _varianceIter;
};

template class CompressFunctor<float>;
template class CompressFunctor<double>;

}}} //end namespace lsst::meas::multifit

#endif
