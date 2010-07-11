// -*- lsst-c++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
#include "lsst/meas/multifit/footprintUtils.h"

namespace multifit = lsst::meas::multifit;
namespace afwDet = lsst::afw::detection;

/**
 * Construct an lsst::afw::detection::Footprint from an
 * lsst::afw::geom::ellipses::Ellipse. 
 *
 * @note Move to lsst::afw::detection
 */
afwDet::Footprint::Ptr multifit::makeFootprint(
    lsst::afw::geom::ellipses::Ellipse const & ellipse
) {
    lsst::afw::detection::Footprint::Ptr fp(new lsst::afw::detection::Footprint());
    lsst::afw::geom::ellipses::Core::RadialFraction rf(ellipse.getCore());

    lsst::afw::geom::BoxD envelope = ellipse.computeEnvelope();
    int const yEnd = static_cast<int>(envelope.getMaxY()) + 1;
    int const xEnd = static_cast<int>(envelope.getMaxX()) + 1;
    lsst::afw::geom::ExtentD dp(lsst::afw::geom::PointD(0) -ellipse.getCenter());
    for (int y = static_cast<int>(envelope.getMinY()); y<yEnd; ++y) {
        int x = static_cast<int>(envelope.getMinX());
        while (rf(lsst::afw::geom::PointD::make(x,y) + dp) > 1.0) {
            if (x >= xEnd) {
                if (++y >= yEnd) 
                    return fp;
                x = static_cast<int>(envelope.getMinX());
            } else {
                ++x;
            }
        }
        int start = x;
        while (rf(lsst::afw::geom::PointD::make(x,y) + dp) <= 1.0 && x < xEnd) 
            ++x;

        fp->addSpan(y, start, x-1);
    }
    fp->normalize();
    return fp;

#if 0

    afwDet::Footprint::Ptr fp(new afwDet::Footprint());

    lsst::afw::geom::ellipses::Axes axes(ellipse.getCore());
    lsst::afw::geom::BoxD envelope = ellipse.computeEnvelope();
    lsst::afw::geom::ExtentD center(ellipse.getCenter());
    
    int const yMin = static_cast<int>(std::floor(envelope.getMinY()));
    int const yMax = static_cast<int>(std::ceil(envelope.getMaxY()));
    int offsetY = static_cast<int>(envelope.getMinY() - center.getY());  
    double sinTheta = std::sin(axes[lsst::afw::geom::ellipses::Axes::THETA]);
    double cosTheta = std::cos(axes[lsst::afw::geom::ellipses::Axes::THETA]);
    double sinSqTheta = sinTheta*sinTheta;
    double cosSqTheta = cosTheta*cosTheta;    

    double major = axes[lsst::afw::geom::ellipses::Axes::A];
    double minor = axes[lsst::afw::geom::ellipses::Axes::B];
    double majorSq = major*major;
    double minorSq = minor*minor;


    //precompute parts of quadratic equation that do not depend on y
    double aTimes2 = 2*(minorSq*cosSqTheta + majorSq*sinSqTheta);
    double bMult = 2*sinTheta*cosTheta*(minorSq-majorSq);
    double cMult = minorSq*sinSqTheta + majorSq*cosSqTheta;
    double cAdd = - majorSq*minorSq;

    //loop over all rows covering the ellipse
    //solve quadratic equation for min/max x as a function of y
    
    for (int y = yMin; y <= yMax; ++y) {
        double b = bMult * offsetY; 
        double c = cMult * offsetY * offsetY + cAdd;

        double discriminant = b*b - 2 * aTimes2 * c;

        //should usually have at least one solution,
        //but because we are taking floor/ceil of doubles to compute our bbox, 
        //it is possible that the top or bottom row  are empty
        if(discriminant < 0)
            continue;

        double sqRootDiscriminant = std::sqrt(discriminant);
        int x0 = static_cast<int>(std::floor( (-b - sqRootDiscriminant) / (aTimes2) + center.getX()));
        int x1 = static_cast<int>(std::ceil( (-b + sqRootDiscriminant) / (aTimes2) + center.getX()));

        fp->addSpan(y, x0, x1);
        ++offsetY;
    }
    return fp;
#endif
}

/**
 * Using an image mask, clip the footprint to fit on the image, and remove
 * from the footprint any masked pixels
 *
 * This invokes a deep copy of the footprint, even if the image is large enough
 * to accomodate the entire footprint, and no pixels are masked
 */
template <typename MaskPixel> 
afwDet::Footprint::Ptr multifit::clipAndMaskFootprint(
    lsst::afw::detection::Footprint const & footprint,
    typename lsst::afw::image::Mask<MaskPixel>::Ptr const & mask,
    MaskPixel bitmask
) {
    if(bitmask == 0)
        bitmask = ~bitmask;

    lsst::afw::image::BBox const maskBBox(
        mask->getXY0(), mask->getWidth(), mask->getHeight()
    );
    int maskX0 = maskBBox.getX0();
    int maskY0 = maskBBox.getY0();
    int maskX1 = maskBBox.getX1();
    int maskY1 = maskBBox.getY1();
    
    afwDet::Footprint::Ptr fixedFootprint(new afwDet::Footprint()); 

    afwDet::Footprint::SpanList const & oldSpanList(
        footprint.getSpans()
    );

    afwDet::Footprint::SpanList::const_iterator spanIter(
        oldSpanList.begin()
    );
    afwDet::Footprint::SpanList::const_iterator const & spanEnd(
        oldSpanList.end()
    );
    

    int x0, x1, y;
    for( ; spanIter != spanEnd; ++spanIter) {
        afwDet::Span const & span(**spanIter);
        y = span.getY();
        x0 = span.getX0();
        x1 = span.getX1();
        if(y < maskY0 || y > maskY1 || x1 < maskX0 || x0 > maskX1) {
            //span is entirely outside the image mask. cannot be used
            continue;
        }

        //clip the span to be within the mask
        if(x0 < maskX0) x0 = maskX0;
        if(x1 > maskX1) x1 = maskX1;

        //lsst::afw::image iterators are always specified with respect to (0,0)
        //regardless what the image::XY0 is set to.        
        typename lsst::afw::image::Mask<MaskPixel>::x_iterator maskIter = 
            mask->x_at(x0 - maskX0, y - maskY0);

        //loop over all span locations, slicing the span at maskedPixels
        for(int x = x0; x <= x1; ++x) {            
            if((*maskIter & bitmask) != 0) {
                //masked pixel found within span
                if (x > x0) {                    
                    //add beginning of spanto the fixedFootprint
                    //the fixed span contains all the unmasked pixels up to,
                    //but not including this masked pixel
                    fixedFootprint->addSpan(y, x0, x - 1);                
                }
                //set the next Span to start after this pixel
                x0 = x + 1;
            }

            //move to next mask pixel
            ++maskIter;
        }
        
        //add last section of span
        if(x0 <= x1) {
            fixedFootprint->addSpan(y, x0, x1);
        }
    }
   
    fixedFootprint->setRegion(maskBBox);
    return fixedFootprint;
}

//FootprintFunctor to compress an afw::Image to two ndarray::Array using a footprint
template <typename ImagePixel, typename MaskPixel=lsst::afw::image::MaskPixel,
    typename VariancePixel=lsst::afw::image::VariancePixel>
class CompressImageFunctor : public lsst::afw::detection::FootprintFunctor<
        lsst::afw::image::MaskedImage<ImagePixel, MaskPixel, VariancePixel> > {
public:
    typedef lsst::afw::image::MaskedImage<ImagePixel, MaskPixel, VariancePixel> MaskedImage;
    typedef lsst::afw::detection::FootprintFunctor<MaskedImage> FootprintFunctor;

    CompressImageFunctor(
        MaskedImage const & src,
        ndarray::Array<multifit::Pixel, 1, 1> const & imageDest,
        ndarray::Array<multifit::Pixel, 1, 1> const & varianceDest
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
       *_imageIter = static_cast<multifit::Pixel>(loc.image());
       *_varianceIter = static_cast<multifit::Pixel>(loc.variance());
       ++_imageIter;
       ++_varianceIter;
    }

private:
    ndarray::Array<multifit::Pixel, 1, 1> const & _imageDest;
    ndarray::Array<multifit::Pixel, 1, 1> const & _varianceDest;
    ndarray::Array<multifit::Pixel, 1, 1>::Iterator _imageIter, _varianceIter;
};

template class CompressImageFunctor<float>;
template class CompressImageFunctor<double>;

/**
 * Compress an lsst::afw::image::MaskedImage into two ndarray::Array
 *
 * This is a deep copy of all pixels of maskedImage that are in footprint
 * imageDest and varianceDest must be initialized externally, prior to calling
 * this function, and must be of sie footprint.getNpix()
 *
 * @param footprint defines all the desired pixels to copy
 * @param maskedImage the source image to compress
 * @param imageDest destination for the image pixel data
 * @param varianceDest destination for the variance pixel data
 */
template <typename ImagePixel, typename MaskPixel, typename VariancePixel>
void multifit::compressImage(
    lsst::afw::detection::Footprint const & footprint,
    lsst::afw::image::MaskedImage<ImagePixel, MaskPixel, VariancePixel> const & maskedImage,
    ndarray::Array<Pixel, 1, 1> const & imageDest,
    ndarray::Array<Pixel, 1, 1> const & varianceDest
) {
    CompressImageFunctor<ImagePixel, MaskPixel, VariancePixel> functor(
        maskedImage, imageDest, varianceDest
    ); 
    functor.apply(footprint);
}


//FootprintFunctor to expand two ndarray::Array to an afw::Image using a footprint
template <typename ImagePixel, typename MaskPixel=lsst::afw::image::MaskPixel, 
    typename VariancePixel=lsst::afw::image::VariancePixel>
class  ExpandImageFunctor : public lsst::afw::detection::FootprintFunctor<
        lsst::afw::image::MaskedImage<ImagePixel, MaskPixel, VariancePixel> > {
public:
    typedef lsst::afw::image::MaskedImage<ImagePixel, MaskPixel, VariancePixel> MaskedImage;
    typedef lsst::afw::detection::FootprintFunctor<MaskedImage> FootprintFunctor;

    ExpandImageFunctor(
        MaskedImage & dest,
        ndarray::Array<multifit::Pixel const, 1, 1> const & imageSrc,
        ndarray::Array<multifit::Pixel const, 1, 1> const & varianceSrc
    ) : FootprintFunctor(dest),
        _imageSrc(imageSrc),
        _varianceSrc(varianceSrc)        
    {}

    virtual void reset(lsst::afw::detection::Footprint const & footprint) {
        if(_imageSrc.getSize<0>() != footprint.getNpix() || 
            _varianceSrc.getSize<0>() != footprint.getNpix()) {
            throw LSST_EXCEPT(lsst::pex::exceptions::LengthErrorException,
                "Source vectors are not correct length");    
        }        
        _imageIter = _imageSrc.begin();
        _varianceIter = _varianceSrc.begin(); 
    }
    virtual void operator()(
        typename MaskedImage::xy_locator loc,
        int x,
        int y
    ) {
       loc.image() = static_cast<ImagePixel>(*_imageIter);
       loc.variance() = static_cast<VariancePixel>(*_varianceIter);
       ++_imageIter;
       ++_varianceIter;
    }

private:
    ndarray::Array<multifit::Pixel const, 1, 1> const & _imageSrc;
    ndarray::Array<multifit::Pixel const, 1, 1> const & _varianceSrc;
    ndarray::Array<multifit::Pixel const, 1, 1>::Iterator _imageIter, _varianceIter; 
};

template class ExpandImageFunctor<float>;
template class ExpandImageFunctor<double>;

/**
 * Expand two ndarray::Array to a full lsst::afw::image::MaskedImage
 *
 * This is a deep copy of all values in both inout arrays, the footprint
 * determines what pixel in the image corresponds to each position in the array
 * 
 * imageSrc and varianceSrc must be initialized externally, prior to calling
 * this function, and must be of sie footprint.getNpix()
 *
 * maskedImage must be initialized externally, and completely cover all pixels
 * in the footprint. Pixels contained in the footprint will be set by this
 * function. Those pixels outisde the footprint will not be modified.
 *
 * @param footprint maps array index to image location
 * @param maskedImage destination of the copy. only those pixels contained in
 * the footpint will be touched.
 * @param imageSrc source for the image pixel data
 * @param varianceSrc source for the variance pixel data
 */
template <typename ImagePixel, typename MaskPixel, typename VariancePixel>
void multifit::expandImage(
    lsst::afw::detection::Footprint const & footprint,
    lsst::afw::image::MaskedImage<ImagePixel, MaskPixel, VariancePixel> & maskedImage,
    ndarray::Array<Pixel const, 1, 1> const & imageSrc,
    ndarray::Array<Pixel const, 1, 1> const & varianceSrc
) {
    ExpandImageFunctor<ImagePixel, MaskPixel, VariancePixel> functor(
        maskedImage, imageSrc, varianceSrc
    ); 
    functor.apply(footprint);
}



//explicit template instantiations
template afwDet::Footprint::Ptr multifit::clipAndMaskFootprint<lsst::afw::image::MaskPixel>(
    afwDet::Footprint const &,
    lsst::afw::image::Mask<lsst::afw::image::MaskPixel>::Ptr const &,
    lsst::afw::image::MaskPixel 
);

template void multifit::compressImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>(
    lsst::afw::detection::Footprint const &,
    lsst::afw::image::MaskedImage<float> const &,
    ndarray::Array<Pixel, 1, 1> const &,
    ndarray::Array<Pixel, 1, 1> const &
);
template void multifit::compressImage<double, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>(
    lsst::afw::detection::Footprint const &,
    lsst::afw::image::MaskedImage<double> const &,
    ndarray::Array<Pixel, 1, 1> const &,
    ndarray::Array<Pixel, 1, 1> const &
);


template void multifit::expandImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>(
    lsst::afw::detection::Footprint const & footprint,
    lsst::afw::image::MaskedImage<float> &,
    ndarray::Array<Pixel const, 1, 1> const & imageSrc,
    ndarray::Array<Pixel const, 1, 1> const & varianceSrc
);
template void multifit::expandImage<double, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>(
    lsst::afw::detection::Footprint const & footprint,
    lsst::afw::image::MaskedImage<double> &,
    ndarray::Array<Pixel const, 1, 1> const & imageSrc,
    ndarray::Array<Pixel const, 1, 1> const & varianceSrc
);



