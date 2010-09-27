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
 
/**
 * @file
 * Implementation of ModelEvaluator
 */
#include "lsst/meas/multifit/ModelEvaluator.h"
#include "lsst/meas/multifit/matrices.h"
#include "lsst/meas/multifit/footprintUtils.h"

#include <iostream>
namespace multifit = lsst::meas::multifit;
namespace afwImg = lsst::afw::image;
namespace afwDet = lsst::afw::detection;
/**
 * Set the list of exposures used to evaluate the model
 *
 * This is an atomic operation which resets the state of this ModelEvaluator 
 * completely. The ModelEvaluator will not be properly initialized until after
 * this function is called.
 *
 * For each exposure in the list, a projection footprint of the model is 
 * computed. If the projection footprint has more than \c getNMinPix pixels
 * which fall within the bounding box of the exposure, then a projection is
 * generated for that exposure.
 *
 * The pixel threshold can be set on construction or by calling setNMinPix
 *
 * Data and variance vectors are constructed by concactenating all the
 * contributing pixels from each projection.
 *
 * @sa getNMinPix
 * @sa setNMinPix
 */
template<typename ExposureT>
void multifit::ModelEvaluator::setExposures(    
    std::list<ExposureT> const & exposureList
) { 
    typedef std::list<ExposureT> ExposureList;
    typedef typename ExposureList::const_iterator ExposureIterator;
    typedef typename ExposureT::MaskedImageT MaskedImageT;

    std::list<MaskedImageT> imageList;
    std::list<lsst::afw::geom::AffineTransform> transformList;
    std::list<CONST_PTR(lsst::afw::detection::Psf)> psfList;

    for(ExposureIterator i(exposureList.begin()), end(exposureList.end()); i != end; ++i) {
        imageList.push_back(i->getMaskedImage());
        lsst::afw::image::Wcs::Ptr wcs = i->getWcs();
        transformList.push_back(wcs->linearizeSkyToPixel(_model->computePosition()));
        psfList.push_back(i->getPsf());
    }

    setData<MaskedImageT>(imageList, psfList, transformList);
}

template <typename MaskedImageT>
void multifit::ModelEvaluator::setData(
    std::list<MaskedImageT> const & imageList,
    std::list<lsst::afw::detection::Psf::ConstPtr> const & psfList,
    std::list<lsst::afw::geom::AffineTransform> const & skyToPixelTransformList
) {
    if(imageList.size() != skyToPixelTransformList.size() || imageList.size() != psfList.size()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LengthErrorException,
            "Input data lists must contain same number of elements"
        );
    }

    typedef typename MaskedImageT::Mask Mask;
    typedef typename Mask::Pixel MaskPixel;

    _projectionList.clear();
    _validProducts = 0;
    
    int nLinear = getLinearParameterSize();
    int nNonlinear = getNonlinearParameterSize();

    int pixSum = 0;
    
    //exposures which contain fewer than _nMinPix pixels will be rejected
    //construct a list containing only those exposure which were not rejected
    std::list<MaskedImageT> goodImageList;

    MaskPixel bitmask = Mask::getPlaneBitMask("BAD") | 
        Mask::getPlaneBitMask("INTRP") | Mask::getPlaneBitMask("SAT") | 
        Mask::getPlaneBitMask("CR") | Mask::getPlaneBitMask("EDGE");

    // loop to create projections
    
    typename std::list<MaskedImageT>::const_iterator iImage(imageList.begin()), endImage(imageList.end());
    std::list<lsst::afw::geom::AffineTransform>::const_iterator iTransform(
        skyToPixelTransformList.begin()
    );
    std::list<CONST_PTR(lsst::afw::detection::Psf)>::const_iterator iPsf(
        psfList.begin()
    );
    for( ; iImage != endImage; ++iImage, ++iTransform, ++iPsf) {
        afwDet::Footprint::Ptr projectionFp = _model->computeProjectionFootprint(*iPsf, *iTransform);
        afwDet::Footprint::Ptr fixedFp = clipAndMaskFootprint<MaskPixel>(
            *projectionFp, 
            *iImage->getMask(),
            bitmask
        );
        //ignore exposures with too few contributing pixels        
        if (fixedFp->getNpix() > _nMinPix) {
            _projectionList.push_back(
                _model->makeProjection(*iPsf, *iTransform, fixedFp)
            );
            pixSum += fixedFp->getNpix();
            goodImageList.push_back(*iImage);
        }
    }

    //  allocate matrix buffers
    ndarray::shallow(_dataVector) = ndarray::allocate<Allocator>(ndarray::makeVector(pixSum));
    ndarray::shallow(_varianceVector) = ndarray::allocate<Allocator>(ndarray::makeVector(pixSum));
    ndarray::shallow(_modelImageBuffer) = ndarray::allocate<Allocator>(ndarray::makeVector(pixSum));
    ndarray::shallow(_linearDerivativeBuffer) = ndarray::allocate<Allocator>(
        ndarray::makeVector(nLinear, pixSum)
    );
    ndarray::shallow(_nonlinearDerivativeBuffer) = ndarray::allocate<Allocator>(
        ndarray::makeVector(nNonlinear, pixSum)
    );    
    
    int nPix;
    int pixelStart = 0, pixelEnd;
    
    iImage = goodImageList.begin();
    //loop to assign matrix buffers to each projection Frame
    for(ProjectionIterator i(_projectionList.begin()), end(_projectionList.end()); 
        i != end; ++i, ++iImage
    ) {
        ModelProjection & projection(**i);
        nPix = projection.getFootprint()->getNpix();
        pixelEnd = pixelStart + nPix;

        // compress the exposure using the footprint
        compressImage(
            *projection.getFootprint(), 
            *iImage, 
            _dataVector[ndarray::view(pixelStart, pixelEnd)], 
            _varianceVector[ndarray::view(pixelStart, pixelEnd)] 
        );

        //set modelImage buffer
        projection.setModelImageBuffer(
            _modelImageBuffer[ndarray::view(pixelStart, pixelEnd)]
        );
        
        //set linear buffer
        projection.setLinearParameterDerivativeBuffer(
            _linearDerivativeBuffer[ndarray::view()(pixelStart, pixelEnd)]
        );
        //set nonlinear buffer
        projection.setNonlinearParameterDerivativeBuffer(
            _nonlinearDerivativeBuffer[ndarray::view()(pixelStart, pixelEnd)]
        );

        pixelStart = pixelEnd;
    }
    

    VectorMap varianceMap (_varianceVector.getData(), getNPixels(), 1);
    _sigma = varianceMap.cwise().sqrt(); 

}

/**
 * Compute the value of the model at every contributing pixel of every exposure
 *
 * @sa ModelProjection::computeModelImage
 */
Eigen::Matrix<multifit::Pixel, Eigen::Dynamic, 1> const & 
multifit::ModelEvaluator::computeModelImage() {
    if(!(_validProducts & MODEL_IMAGE)) {
        for( ProjectionIterator i(_projectionList.begin()), end(_projectionList.end());
             i  != end; ++i
        ) {
            (*i)->computeModelImage();
        }
        VectorMap map(_modelImageBuffer.getData(), getNPixels(),1);
        _modelImage = map.cwise()/_sigma;
        _validProducts |= MODEL_IMAGE;
    }    
    return _modelImage;
}

/**
 * Compute the derivative of the model with respect to its linear parameters
 *
 * @sa ModelProjection::computeLinearParameterDerivative
 */
Eigen::Matrix<multifit::Pixel, Eigen::Dynamic, Eigen::Dynamic> const & 
multifit::ModelEvaluator::computeLinearParameterDerivative() {
    if (!(_validProducts & LINEAR_PARAMETER_DERIVATIVE)) {
        for( ProjectionIterator i(_projectionList.begin()), 
             end(_projectionList.end()); i  != end; ++i
        ) {
            (*i)->computeLinearParameterDerivative();
        }
        MatrixMap map(_linearDerivativeBuffer.getData(), getNPixels(), getLinearParameterSize());
        _linearDerivative = map;
        for(int i =0; i < getLinearParameterSize(); ++i) {
            _linearDerivative.col(i).cwise() /= _sigma;
        }
        _validProducts |= LINEAR_PARAMETER_DERIVATIVE;
    }    
    return _linearDerivative;
}

/**
 * Compute the derivative of the model with respect to its nonlinear parameters
 *
 * @sa ModelProjection::computeNonlinearParameterDerivative
 */
Eigen::Matrix<multifit::Pixel, Eigen::Dynamic, Eigen::Dynamic> const & 
multifit::ModelEvaluator::computeNonlinearParameterDerivative() {
    if(!(_validProducts & NONLINEAR_PARAMETER_DERIVATIVE)) {
        for( ProjectionIterator i(_projectionList.begin()), 
             end(_projectionList.end()); i  != end; ++i
        ) {
            (*i)->computeNonlinearParameterDerivative();
        }
        MatrixMap map(_nonlinearDerivativeBuffer.getData(), getNPixels(), getNonlinearParameterSize());
        _nonlinearDerivative = map;
        for(int i =0; i < getLinearParameterSize(); ++i) {
            _nonlinearDerivative.col(i).cwise() /= _sigma;
        }
        _validProducts |= NONLINEAR_PARAMETER_DERIVATIVE;
    }    
    return _nonlinearDerivative;
}


template void multifit::ModelEvaluator::setExposures<
        afwImg::Exposure<float, afwImg::MaskPixel, afwImg::VariancePixel>
> (
    std::list<afwImg::Exposure<float> > const &
);
template void multifit::ModelEvaluator::setExposures<
        afwImg::Exposure<double, afwImg::MaskPixel, afwImg::VariancePixel>
> (
    std::list<afwImg::Exposure<double> > const &
);

template void multifit::ModelEvaluator::setData<
        afwImg::MaskedImage<double, afwImg::MaskPixel, afwImg::VariancePixel>
> (
    std::list<afwImg::MaskedImage<double> > const &,
    std::list<CONST_PTR(lsst::afw::detection::Psf)> const &,
    std::list<lsst::afw::geom::AffineTransform> const &
);
template void multifit::ModelEvaluator::setData<
        afwImg::MaskedImage<float, afwImg::MaskPixel, afwImg::VariancePixel>
> (
    std::list<afwImg::MaskedImage<float> > const &,
    std::list<CONST_PTR(lsst::afw::detection::Psf)> const &,
    std::list<lsst::afw::geom::AffineTransform> const &
);
