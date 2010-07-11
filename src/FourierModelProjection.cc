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
 
#include "lsst/meas/multifit/FourierModelProjection.h"
#include <ndarray/fft.hpp>
#include <iostream>

#include "lsst/afw/image/Utils.h"
#include "lsst/afw/geom.h"


namespace multifit = lsst::meas::multifit;

/**
 * Manager of intermdiate products for doing shift operations
 *
 * @note Internal implementation detail of FourierModelProjection
 */
class multifit::FourierModelProjection::Shifter : private boost::noncopyable {
public:
    void handleNonlinearParameterChange() { _valid = false; }

    void apply(
        ndarray::FourierArray<Pixel,3,3>::Iterator iter, 
        ndarray::FourierArray<Pixel,3,3>::Iterator const & end
    ) {
        if (!_valid) {
            lsst::afw::geom::Point2D translation(
                _parent->_getPsfPosition() - lsst::afw::geom::PointD(_parent->_outerBBox.getMin())
            );

            _factor = 1.0;

            ndarray::shift(                
                ndarray::makeVector(translation.getY(), translation.getX()), 
                _factor
            );

            _valid = true;
        }

        for (; iter != end; ++iter) {
            *iter *= _factor;
        }
    }

    explicit Shifter(FourierModelProjection * parent) :
        _parent(parent), _valid(false),
        _factor(FFT::initializeK(
            ndarray::makeVector(
                _parent->_outerBBox.getHeight(), 
                _parent->_outerBBox.getWidth()
            )
        ))
    { }

private:
    FourierModelProjection * _parent;
    bool _valid;
    ndarray::FourierArray<Pixel,2,2> _factor;
};

/**
 * Manager of intermdiate products for computing the derivative with respect to
 * linear parameters
 *
 * @note Internal implementation detail of FourierModelProjection
 */
class multifit::FourierModelProjection::LinearMatrixHandler : boost::noncopyable {
public:
    void handleLinearParameterChange() { _imageValid = false; }
    void handleNonlinearParameterChange() { 
        _matrixValid = false; _imageValid = false; 
    }
    
    ndarray::Array<Pixel const,3,1> computeLinearParameterDerivative() {
        if (!_matrixValid) {
            _kLPD = _parent->_getMorphologyProjection()->computeLinearParameterDerivative();
            if(!_parent->_shifter)
                _parent->_shifter.reset(new Shifter(_parent));
            _parent->_shifter->apply(_kLPD.begin(),_kLPD.end());
            _unconvolvedMatrix = _kLPD;
            _parent->_applyKernel(_kLPD.begin(),_kLPD.end());
            _ifft->execute();
            _matrixValid = true;

        }

        return _finalLPD;
    }

    ndarray::FourierArray<Pixel const,2,2> computeUnconvolvedImage() {
        computeLinearParameterDerivative();
        if (!_imageValid) {
            getCompressedVectorView(_unconvolvedImage.getBase()) = 
                getCompressedMatrixView(_unconvolvedMatrix.getBase()) * 
                (_parent->getModel()->getLinearParameters());
            _imageValid = true;
        }

        return _unconvolvedImage;
    }

    explicit LinearMatrixHandler(FourierModelProjection * parent) 
      : _matrixValid(false), 
        _imageValid(false), 
        _parent(parent)
    {
        _ifft = FFT::planMultiplexInverse(
            ndarray::makeVector(
                _parent->getLinearParameterSize(),
                _parent->_outerBBox.getHeight(),
                _parent->_outerBBox.getWidth()
            ),
            _kLPD, 
            _xLPD
        );
        ndarray::shallow(_unconvolvedMatrix) = FFT::initializeK(
            _xLPD.getShape()
        );
        ndarray::shallow(_unconvolvedImage) = FFT::initializeK(
            ndarray::makeVector(
                _parent->_outerBBox.getHeight(), 
                _parent->_outerBBox.getWidth()
            )
        );

        ndarray::shallow(_finalLPD) = window(_xLPD, _parent->_innerBBox);
    }

private:
    bool _matrixValid;
    bool _imageValid;
    FourierModelProjection * _parent;
    FFT::Ptr _ifft;
    ndarray::FourierArray<Pixel,3,3> _unconvolvedMatrix;
    ndarray::FourierArray<Pixel,2,2> _unconvolvedImage;
    ndarray::FourierArray<Pixel,3,3> _kLPD;
    ndarray::Array<Pixel,3,3> _xLPD;
    ndarray::Array<Pixel,3,1> _finalLPD;
};

/**
 * Manager of intermdiate products for computing the derivative with respect to
 * nonlinear parameters
 *
 * @note Internal implementation detail of FourierModelProjection
 */
class multifit::FourierModelProjection::NonlinearMatrixHandler : boost::noncopyable {
public:
    void handleParameterChange() { _valid = false; }
    
    ndarray::Array<Pixel const,3,1> computeTranslationDerivative() {
        recompute();
        return _finalTD; 
    }
    ndarray::Array<Pixel const,3,1> computeProjectedParameterDerivative() {
        recompute();
        return _finalPPD; 
    }

    explicit NonlinearMatrixHandler(FourierModelProjection * parent) :
        _valid(false), _parent(parent)
    {
        _ifft = FFT::planMultiplexInverse(
            ndarray::makeVector(
                _parent->getNonlinearParameterSize(),
                _parent->_outerBBox.getHeight(),
                _parent->_outerBBox.getWidth()
            ),
            _kFull, 
            _xFull
        );
        ndarray::Array<Pixel,3,1> finalFull(window(_xFull, _parent->_innerBBox));
        if (_parent->hasTranslationDerivative()) {
            ndarray::shallow(_kTD) = _kFull[ndarray::view(0, 2)];
            ndarray::shallow(_finalTD) = finalFull[ndarray::view(0, 2)];
        }
        if (_parent->hasProjectedParameterDerivative()) {
            ndarray::shallow(_kPPD) = _kFull[ndarray::view(2, _kFull.getSize<0>())];
            ndarray::shallow(_finalPPD) = finalFull[ndarray::view(2, finalFull.getSize<0>())];
        }
    }

private:
    void recompute() {
        if (_valid) 
            return;
        if (_parent->hasTranslationDerivative()) {
            _kTD[0] = _kTD[1] = _parent->_linearMatrixHandler->computeUnconvolvedImage();
            _parent->_applyKernel(_kTD.begin(),_kTD.end());
            ndarray::differentiate(1, _kTD[0]);
            ndarray::differentiate(0, _kTD[1]);
        }
        if (_parent->hasProjectedParameterDerivative()) {
            _kPPD = _parent->_getMorphologyProjection()->computeProjectedParameterDerivative(); 
            if(!_parent->_shifter)
                _parent->_shifter.reset(new Shifter(_parent));
            _parent->_shifter->apply(_kPPD.begin(), _kPPD.end());
            _parent->_applyKernel(_kPPD.begin(), _kPPD.end());
        }
        _ifft->execute();
        _valid = true;
    }

    bool _valid;
    FourierModelProjection * _parent;
    FFT::Ptr _ifft;
    ndarray::FourierArray<Pixel,3,3> _kFull;
    ndarray::FourierArray<Pixel,3,3> _kTD;  // translation derivative
    ndarray::FourierArray<Pixel,3,3> _kPPD; // projected parameter derivative
    ndarray::Array<Pixel,3,3> _xFull;
    ndarray::Array<Pixel,3,1> _finalTD;
    ndarray::Array<Pixel,3,1> _finalPPD;
};

/**
 * Manager of intermdiate products for computing the derivative with respect to
 * psf parameters
 *
 * @note Internal implementation detail of FourierModelProjection
 */
class multifit::FourierModelProjection::PsfMatrixHandler : boost::noncopyable {
public:
    void handleParameterChange() { _valid = false; }
    
    ndarray::Array<Pixel const,3,1> computePsfParameterDerivative() {
        if (_valid) 
            return _final;
        
        lsst::afw::math::FourierCutout::Ptr fourierCutout = 
            _parent->_localKernel->getFourierDerivatives()[0];

        ndarray::Array<std::complex<Pixel>,3,3> externalImg(
            ndarray::external(
                fourierCutout->begin(),
                ndarray::makeVector(
                    _parent->_localKernel->getNParameters(), 
                    fourierCutout->getFourierHeight(),
                    fourierCutout->getFourierWidth()
                ),
                ndarray::makeVector(
                    fourierCutout->getFourierSize(),
                    fourierCutout->getFourierWidth(),
                    1
                ),
                fourierCutout->getOwner()
            )
        );
        ndarray::FourierArray<Pixel, 3, 3> fourierDerivative(
            fourierCutout->getImageWidth(),
            externalImg
        );
        _k = fourierDerivative;
        _parent->_applyKernel(
            _k.begin(), _k.end(),
            _parent->_linearMatrixHandler->computeUnconvolvedImage()
        );
            
        _ifft->execute();
        _valid = true;
        return _final;
    }

    explicit PsfMatrixHandler(FourierModelProjection * parent) :
        _valid(false), _parent(parent) 
    {
        _ifft = FFT::planMultiplexInverse(
            ndarray::makeVector(
                _parent->getPsfParameterSize(),
                _parent->_outerBBox.getHeight(),
                _parent->_outerBBox.getWidth()
            ),
            _k, _x
        );
        ndarray::shallow(_final) = window(_x, _parent->_innerBBox);
    }

private:
    bool _valid;
    FourierModelProjection * _parent;
    FFT::Ptr _ifft;
    ndarray::FourierArray<Pixel,3,3> _k;
    ndarray::Array<Pixel,3,3> _x;
    ndarray::Array<Pixel,3,1> _final;
};

// ----------------------------- FourierModelProjection ------------------------


int const multifit::FourierModelProjection::getPsfParameterSize() const {
    if(!_localKernel)
        return 0;
    return _localKernel->getNParameters();
}

/**
 * Compute a locally evaluated ConvolutionVisitor from the PSF's
 * kernel, and apply it to the model. This operation invalidates all
 * intermediadiary products. 
 */
void multifit::FourierModelProjection::_convolve(
    PsfConstPtr const & psf
) { 
    if(!psf || !psf->getKernel()) 
        return;

    _localKernel = psf->getKernel()->computeFourierLocalKernel(
        _getPsfPosition()
    );
    
    if(_psfMatrixHandler){
        _psfMatrixHandler.reset(new PsfMatrixHandler(this));
    }
    if(_linearMatrixHandler) {
        _linearMatrixHandler->handleNonlinearParameterChange();
    }
    if (_nonlinearMatrixHandler) {
        _nonlinearMatrixHandler->handleParameterChange();
    }
}

void multifit::FourierModelProjection::_computeLinearParameterDerivative(
    ndarray::Array<Pixel,2,1> const & output
) { 
    if(!_linearMatrixHandler)
        _linearMatrixHandler.reset(new LinearMatrixHandler(this));
    
    _wf->compress(
        _linearMatrixHandler->computeLinearParameterDerivative(),
        output
    );
}

void multifit::FourierModelProjection::_computePsfParameterDerivative(
    ndarray::Array<Pixel,2,1> const & output
) {
    if(!_psfMatrixHandler)
        _psfMatrixHandler.reset(new PsfMatrixHandler(this));

    _wf->compress(
        _psfMatrixHandler->computePsfParameterDerivative(),
        output
    );
}

void multifit::FourierModelProjection::_computeTranslationDerivative(
    ndarray::Array<Pixel,2,1> const & output
) {
    if(!_nonlinearMatrixHandler)
        _nonlinearMatrixHandler.reset(new NonlinearMatrixHandler(this));
    
    _wf->compress(
        _nonlinearMatrixHandler->computeTranslationDerivative(),
        output
    );
}

void multifit::FourierModelProjection::_computeProjectedParameterDerivative(
    ndarray::Array<Pixel,2,1> const & output
) {
    if(!_nonlinearMatrixHandler)
        _nonlinearMatrixHandler.reset(new NonlinearMatrixHandler(this));
    _wf->compress(
        _nonlinearMatrixHandler->computeProjectedParameterDerivative(), 
        output
    );
}

void multifit::FourierModelProjection::_handleLinearParameterChange() {
    ComponentModelProjection::_handleLinearParameterChange();

    if (_getMorphologyProjection()->getDimensions() != _outerBBox.getDimensions()) {
        _setDimensions();
    } else {
        _linearMatrixHandler->handleLinearParameterChange();
        if (_nonlinearMatrixHandler) _nonlinearMatrixHandler->handleParameterChange();
        if (_psfMatrixHandler) _psfMatrixHandler->handleParameterChange();
    }
}

void multifit::FourierModelProjection::_handleNonlinearParameterChange() {
    ComponentModelProjection::_handleNonlinearParameterChange();
    if (_getMorphologyProjection()->getDimensions() != _outerBBox.getDimensions()) {
        _setDimensions();
    } else {
        _linearMatrixHandler->handleNonlinearParameterChange();
        if (_nonlinearMatrixHandler) _nonlinearMatrixHandler->handleParameterChange();
        if (_psfMatrixHandler) _psfMatrixHandler->handleParameterChange();
    }
    if(!_shifter)
        _shifter.reset(new Shifter(this));

    _shifter->handleNonlinearParameterChange();
}

multifit::FourierModelProjection::FourierModelProjection(
    ComponentModel::ConstPtr const & model,
    PsfConstPtr const & psf,
    WcsConstPtr const & wcs,
    FootprintConstPtr const & footprint
) : ComponentModelProjection(model,psf,wcs,footprint),
    _localKernel(), _wf(), 
    _outerBBox(lsst::afw::geom::Point2I(), lsst::afw::geom::Extent2I()),
    _innerBBox(lsst::afw::geom::Point2I(), lsst::afw::geom::Extent2I())
{
    _convolve(psf);
    _setDimensions();
}

/**
 * Determine size of all arrays
 *
 * Called at initialization, and (rarely) when parameters change significantly
 * 
 */
void multifit::FourierModelProjection::_setDimensions() {
    lsst::afw::geom::Extent2I dimensions = getMorphologyProjection()->getDimensions();
    lsst::afw::geom::PointD centerOnExposure = _getPsfPosition();

    lsst::afw::geom::Point2I bboxMin = lsst::afw::geom::Point2I::make(
        int(std::floor(centerOnExposure.getX() - dimensions.getX()/2)),
        int(std::floor(centerOnExposure.getY() - dimensions.getY()/2))
    );

    _outerBBox = lsst::afw::geom::BoxI(bboxMin, dimensions);
    
    // At this point, _innerBBox is defined relative to exposure    
    lsst::afw::geom::Extent2I padding = getMorphologyProjection()->getPadding();
    _innerBBox = _outerBBox;
    _innerBBox.grow(-padding);
    
    lsst::afw::geom::BoxI footprintBBox = lsst::afw::geom::convertToGeom(getFootprint()->getBBox());
    _innerBBox.clip(footprintBBox);
    if(_innerBBox.isEmpty()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "FourierModelProjection's footprint bounding box, and computed bounding box do not overlap"
        );
    }
    _wf = boost::make_shared<WindowedFootprint>(*getFootprint(), _innerBBox);

    // Shift _innerBBox to be defined relative to _outerBBox.
    _innerBBox.shift(lsst::afw::geom::Point2I() - _outerBBox.getMin());

    _localKernel->setDimensions(_outerBBox.getWidth(), _outerBBox.getHeight());
    _shifter.reset(new Shifter(this));
    _linearMatrixHandler.reset(new LinearMatrixHandler(this));
    if (_nonlinearMatrixHandler)
        _nonlinearMatrixHandler.reset(new NonlinearMatrixHandler(this));
    if (_psfMatrixHandler)
        _psfMatrixHandler.reset(new PsfMatrixHandler(this));
}

/**
 * Utility function for convolving a collection of fourier-space images by a
 * fourier-space kernel image
 */
void multifit::FourierModelProjection::_applyKernel(
    ndarray::FourierArray<Pixel,3,3>::Iterator iter, 
    ndarray::FourierArray<Pixel,3,3>::Iterator const & end,
    ndarray::FourierArray<Pixel const,2,2> const & kernel
) const {

    //iterate over collection of fourier-space images
    for (; iter != end; ++iter) {
        //convolved by kernel-image
        *iter *= kernel;
    }    
}

/**
 * Utility function for apply the projection's PSF's kernel-image to a
 * collection of fourier-space images
 */
void multifit::FourierModelProjection::_applyKernel(
    ndarray::FourierArray<Pixel,3,3>::Iterator iter, 
    ndarray::FourierArray<Pixel,3,3>::Iterator const & end
) const {

    lsst::afw::math::FourierCutout::Ptr fourierCutout = 
        _localKernel->getFourierImage();
  
    //Create an Array over the fourierCutout allocated image
    ndarray::Array<std::complex<Pixel>, 2, 2> externalImg(
        ndarray::external(
            fourierCutout->begin(),
            ndarray::makeVector(
                fourierCutout->getFourierHeight(),
                fourierCutout->getFourierWidth()
            ),
            ndarray::makeVector(fourierCutout->getFourierWidth(), 1),
            fourierCutout->getOwner()
        )
    );
    //Create a FourierArray using the above Array
    ndarray::FourierArray<Pixel, 2, 2> kernelImage(
        fourierCutout->getImageWidth(),
        externalImg
    );


    //delegate to generic overload of this function
    _applyKernel(iter, end, kernelImage);
}

multifit::FourierModelProjection::~FourierModelProjection() {}
