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
namespace afwImg = lsst::afw::image;
namespace afwGeom = lsst::afw::geom;

/**
 * Construct an afwDet::Footprint from an
 * afwGeom::ellipses::Ellipse. 
 *
 * @note Move to afwDet
 */
afwDet::Footprint::Ptr multifit::makeFootprint(
    afwGeom::ellipses::Ellipse const & ellipse
) {
    afwDet::Footprint::Ptr fp(new afwDet::Footprint());
    afwGeom::ellipses::Core::RadialFraction rf(ellipse.getCore());

    afwGeom::BoxD envelope = ellipse.computeEnvelope();
    int const yEnd = static_cast<int>(envelope.getMaxY()) + 1;
    int const xEnd = static_cast<int>(envelope.getMaxX()) + 1;
    afwGeom::ExtentD dp(afwGeom::PointD(0) -ellipse.getCenter());
    for (int y = static_cast<int>(envelope.getMinY()); y<yEnd; ++y) {
        int x = static_cast<int>(envelope.getMinX());
        while (rf(afwGeom::PointD::make(x,y) + dp) > 1.0) {
            if (x >= xEnd) {
                if (++y >= yEnd) 
                    return fp;
                x = static_cast<int>(envelope.getMinX());
            } else {
                ++x;
            }
        }
        int start = x;
        while (rf(afwGeom::PointD::make(x,y) + dp) <= 1.0 && x < xEnd) 
            ++x;

        fp->addSpan(y, start, x-1);
    }
    fp->normalize();
    return fp;

#if 0

    afwDet::Footprint::Ptr fp(new afwDet::Footprint());

    afwGeom::ellipses::Axes axes(ellipse.getCore());
    afwGeom::BoxD envelope = ellipse.computeEnvelope();
    afwGeom::ExtentD center(ellipse.getCenter());
    
    int const yMin = static_cast<int>(std::floor(envelope.getMinY()));
    int const yMax = static_cast<int>(std::ceil(envelope.getMaxY()));
    int offsetY = static_cast<int>(envelope.getMinY() - center.getY());  
    double sinTheta = std::sin(axes[afwGeom::ellipses::Axes::THETA]);
    double cosTheta = std::cos(axes[afwGeom::ellipses::Axes::THETA]);
    double sinSqTheta = sinTheta*sinTheta;
    double cosSqTheta = cosTheta*cosTheta;    

    double major = axes[afwGeom::ellipses::Axes::A];
    double minor = axes[afwGeom::ellipses::Axes::B];
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
 * Produces a deep copy of the footprint, even if the image is large enough
 * to accomodate the entire footprint, and no pixels are masked
 *
 * @param footprint original footprint
 * @param mask used to determine bbox of valid region, and remove masked pixels
 * @param bitmask specifies which mask planes are relevant. By default all mask
 *      planes are considered relevant (~0x0). To ignore masking, specify 0
 */
template <typename MaskPixel> 
afwDet::Footprint::Ptr multifit::clipAndMaskFootprint(
    afwDet::Footprint const & footprint,
    afwImg::Mask<MaskPixel> const & mask,
    MaskPixel const & bitmask
) {
    afwImg::BBox maskBBox(
        mask.getXY0(), mask.getWidth(), mask.getHeight()
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

        //afwImg iterators are always specified with respect to (0,0)
        //regardless what the image::XY0 is set to.        
        typename afwImg::Mask<MaskPixel>::const_x_iterator maskIter = 
            mask.x_at(x0 - maskX0, y - maskY0);

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

namespace {

//FootprintFunctor to compress an afw::MaskedImage to two ndarray::Array using a footprint
template <typename ImageT, typename MaskT=afwImg::MaskPixel, typename VarianceT=afwImg::VariancePixel>
class CompressMaskedImageFunctor : 
    public afwDet::FootprintFunctor<afwImg::MaskedImage<ImageT, MaskT, VarianceT> > {
public:
    typedef afwImg::MaskedImage<ImageT, MaskT, VarianceT> MaskedImage;
    typedef afwDet::FootprintFunctor<MaskedImage> FootprintFunctor;

    CompressMaskedImageFunctor(
        MaskedImage const & src,
        ndarray::Array<multifit::Pixel, 1, 1> const & imageDest,
        ndarray::Array<multifit::Pixel, 1, 1> const & varianceDest
    ) : FootprintFunctor(src),
        _imageDest(imageDest),
        _varianceDest(varianceDest)        
    {}

    virtual void reset(afwDet::Footprint const & footprint) {
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

//FootprintFunctor to compress an afw::Image to an ndarray::Array using a footprint
template <typename ImageT>
class CompressImageFunctor : 
    public afwDet::FootprintFunctor<afwImg::Image<ImageT> > {
public:
    typedef afwImg::Image<ImageT> Image;
    typedef afwDet::FootprintFunctor<Image> FootprintFunctor;

    CompressImageFunctor(
        Image const & src,
        ndarray::Array<multifit::Pixel, 1, 1> const & imageDest
    ) : FootprintFunctor(src),
        _imageDest(imageDest)
    {}

    virtual void reset(afwDet::Footprint const & footprint) {
        if(_imageDest.getSize<0>() != footprint.getNpix()) {
            throw LSST_EXCEPT(lsst::pex::exceptions::LengthErrorException,
                              "Destination vector is not the correct length");    
        }
        _imageIter = _imageDest.begin();
    }
    virtual void operator()(
        typename Image::xy_locator loc,
        int x,
        int y
    ) {
       *_imageIter = static_cast<multifit::Pixel>(*loc);
       ++_imageIter;
    }

private:
    ndarray::Array<multifit::Pixel, 1, 1> const & _imageDest;
    ndarray::Array<multifit::Pixel, 1, 1>::Iterator _imageIter;
};

template class CompressImageFunctor<float>;
template class CompressImageFunctor<double>;

template class CompressMaskedImageFunctor<float>;
template class CompressMaskedImageFunctor<double>;

} //end annonymous namespace

/**
 * Compress an afwImg::MaskedImage into two ndarray::Array
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
template <typename ImageT, typename MaskT, typename VarianceT>
void multifit::compressImage(
    afwDet::Footprint const & footprint,
    afwImg::MaskedImage<ImageT, MaskT, VarianceT> const & maskedImage,
    ndarray::Array<Pixel, 1, 1> const & imageDest,
    ndarray::Array<Pixel, 1, 1> const & varianceDest
) {
    ::CompressMaskedImageFunctor<ImageT, MaskT, VarianceT> functor(
        maskedImage, imageDest, varianceDest
    ); 
    functor.apply(footprint);
}


/**
 * Compress an afwImg::Image into a single ndarray::Array
 *
 * This is a deep copy of all pixels of image that are in footprint.
 * imageDest must be initialized externally, prior to calling
 * this function, and must be of size footprint.getNpix()
 *
 * @param footprint defines all the desired pixels to copy
 * @param image the source image to compress
 * @param imageDest destination for the image pixel data
 */
template <typename ImageT>
void multifit::compressImage(
    afwDet::Footprint const & footprint,
    afwImg::Image<ImageT> const & image,
    ndarray::Array<Pixel, 1, 1> const & imageDest
) {
    ::CompressImageFunctor<ImageT> functor(image, imageDest); 
    functor.apply(footprint);
}


namespace {

//FootprintFunctor to expand two ndarray::Array to an afw::MaskedImage using a footprint
template <typename ImageT, typename MaskT=afwImg::MaskPixel, typename VarianceT=afwImg::VariancePixel>
class  ExpandMaskedImageFunctor : 
    public afwDet::FootprintFunctor<afwImg::MaskedImage<ImageT, MaskT, VarianceT> > {
public:
    typedef afwImg::MaskedImage<ImageT, MaskT, VarianceT> MaskedImage;
    typedef afwDet::FootprintFunctor<MaskedImage> FootprintFunctor;

    ExpandMaskedImageFunctor(
        MaskedImage & dest,
        ndarray::Array<multifit::Pixel const, 1, 1> const & imageSrc,
        ndarray::Array<multifit::Pixel const, 1, 1> const & varianceSrc
    ) : FootprintFunctor(dest),
        _imageSrc(imageSrc),
        _varianceSrc(varianceSrc)        
    {}

    virtual void reset(afwDet::Footprint const & footprint) {
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
       loc.image() = static_cast<ImageT>(*_imageIter);
       loc.variance() = static_cast<VarianceT>(*_varianceIter);
       ++_imageIter;
       ++_varianceIter;
    }

private:
    ndarray::Array<multifit::Pixel const, 1, 1> const & _imageSrc;
    ndarray::Array<multifit::Pixel const, 1, 1> const & _varianceSrc;
    ndarray::Array<multifit::Pixel const, 1, 1>::Iterator _imageIter, _varianceIter; 
};

//FootprintFunctor to expand two ndarray::Array to an afw::Image using a footprint
template <typename ImageT>
class  ExpandImageFunctor : 
    public afwDet::FootprintFunctor<afwImg::Image<ImageT> > {
public:
    typedef afwImg::Image<ImageT> Image;
    typedef afwDet::FootprintFunctor<Image> FootprintFunctor;

    ExpandImageFunctor(
        Image & dest,
        ndarray::Array<multifit::Pixel const, 1, 1> const & imageSrc
    ) : FootprintFunctor(dest),
        _imageSrc(imageSrc)
    {}

    virtual void reset(afwDet::Footprint const & footprint) {
        if(_imageSrc.getSize<0>() != footprint.getNpix()) {
            throw LSST_EXCEPT(lsst::pex::exceptions::LengthErrorException,
                "Source vectors are not correct length");    
        }        
        _imageIter = _imageSrc.begin();
    }
    virtual void operator()(
        typename Image::xy_locator loc,
        int x,
        int y
    ) {
       *loc = static_cast<ImageT>(*_imageIter);
       ++_imageIter;
    }

private:
    ndarray::Array<multifit::Pixel const, 1, 1> const & _imageSrc;
    ndarray::Array<multifit::Pixel const, 1, 1>::Iterator _imageIter; 
};

template class ExpandMaskedImageFunctor<float>;
template class ExpandMaskedImageFunctor<double>;

template class ExpandImageFunctor<float>;
template class ExpandImageFunctor<double>;

} //end annonymous namespace

/**
 * Expand two ndarray::Array to a full afwImg::MaskedImage
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
template <typename ImageT, typename MaskT, typename VarianceT>
void multifit::expandImage(
    afwDet::Footprint const & footprint,
    afwImg::MaskedImage<ImageT, MaskT, VarianceT> & maskedImage,
    ndarray::Array<Pixel const, 1, 1> const & imageSrc,
    ndarray::Array<Pixel const, 1, 1> const & varianceSrc
) {
    ::ExpandMaskedImageFunctor<ImageT, MaskT, VarianceT> functor(
        maskedImage, imageSrc, varianceSrc
    ); 
    functor.apply(footprint);
}

/**
 * Expand an ndarray::Array to a full afwImg::Image
 *
 * This is a deep copy of all values in both inout arrays, the footprint
 * determines what pixel in the image corresponds to each position in the array
 * 
 * imageSrc must be initialized externally, prior to calling
 * this function, and must be of size footprint.getNpix()
 *
 * image must be initialized externally, and completely cover all pixels
 * in the footprint. Pixels contained in the footprint will be set by this
 * function. Those pixels outisde the footprint will not be modified.
 *
 * @param footprint maps array index to image location
 * @param image destination of the copy. only those pixels contained in
 * the footprint will be touched.
 * @param imageSrc source for the image pixel data
 */
template <typename ImageT>
void multifit::expandImage(
    afwDet::Footprint const & footprint,
    afwImg::Image<ImageT> & image,
    ndarray::Array<Pixel const, 1, 1> const & imageSrc
) {
    ::ExpandImageFunctor<ImageT> functor(image, imageSrc); 
    functor.apply(footprint);
}



//explicit template instantiations
template afwDet::Footprint::Ptr multifit::clipAndMaskFootprint<afwImg::MaskPixel>(
    afwDet::Footprint const &,
    afwImg::Mask<afwImg::MaskPixel> const &,
    afwImg::MaskPixel const & 
);

template void multifit::compressImage<float>(
    afwDet::Footprint const &,
    afwImg::MaskedImage<float> const &,
    ndarray::Array<Pixel, 1, 1> const &,
    ndarray::Array<Pixel, 1, 1> const &
);
template void multifit::compressImage<double>(
    afwDet::Footprint const &,
    afwImg::MaskedImage<double> const &,
    ndarray::Array<Pixel, 1, 1> const &,
    ndarray::Array<Pixel, 1, 1> const &
);

template void multifit::expandImage<float>(
    afwDet::Footprint const & footprint,
    afwImg::MaskedImage<float> &,
    ndarray::Array<Pixel const, 1, 1> const & imageSrc,
    ndarray::Array<Pixel const, 1, 1> const & varianceSrc
);
template void multifit::expandImage<double>(
    afwDet::Footprint const & footprint,
    afwImg::MaskedImage<double> &,
    ndarray::Array<Pixel const, 1, 1> const & imageSrc,
    ndarray::Array<Pixel const, 1, 1> const & varianceSrc
);

template void multifit::compressImage<float>(
    afwDet::Footprint const &,
    afwImg::Image<float> const &,
    ndarray::Array<Pixel, 1, 1> const &
);
template void multifit::compressImage<double>(
    afwDet::Footprint const &,
    afwImg::Image<double> const &,
    ndarray::Array<Pixel, 1, 1> const &
);

template void multifit::expandImage<float>(
    afwDet::Footprint const & footprint,
    afwImg::Image<float> &,
    ndarray::Array<Pixel const, 1, 1> const & imageSrc
);
template void multifit::expandImage<double>(
    afwDet::Footprint const & footprint,
    afwImg::Image<double> &,
    ndarray::Array<Pixel const, 1, 1> const & imageSrc
);



