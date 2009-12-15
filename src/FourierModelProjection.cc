#include "lsst/meas/multifit/FourierModelProjection.h"
#include <ndarray/fft.hpp>

#include "lsst/afw/image/Utils.h"
#include "lsst/afw/geom.h"


namespace multifit = lsst::meas::multifit;

// -- multifit::FourierModelProjection::Shifter --------------------------------

class multifit::FourierModelProjection::Shifter : private boost::noncopyable {
public:
    void handleNonlinearParameterChange() { _valid = false; }

    void apply(
        ndarray::FourierArray<Pixel,3,3>::Iterator iter, 
        ndarray::FourierArray<Pixel,3,3>::Iterator const & end
    ) {
        if (!_valid) {
            lsst::afw::geom::Point2D translation(
                _parent->_outerBBox.getMin()
            );
            translation = lsst::afw::geom::Point2D(
                _parent->_getPsfPosition() - translation
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
    {}

private:
    FourierModelProjection * _parent;
    bool _valid;
    ndarray::FourierArray<Pixel,2,2> _factor;
};

// -- multifit::FourierModelProjection::LinearMatrixHandler --------------------

class multifit::FourierModelProjection::LinearMatrixHandler : boost::noncopyable {
public:
    void handleLinearParameterChange() { _imageValid = false; }
    void handleNonlinearParameterChange() { 
        _matrixValid = false; _imageValid = false; 
    }
    
    ndarray::Array<Pixel const,3,1> computeLinearParameterDerivative() {
        if (!_matrixValid) {
            _kLPD = _parent->_getMorphologyProjection()->computeLinearParameterDerivative();
            _kLPD *= _parent->getPhotFactor();
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
                _parent->getModel()->getLinearParameterVector();
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

// -- multifit::FourierModelProjection::NonlinearMatrixHandler -----------------

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
        ndarray::shallow(_kTD) = _kFull[ndarray::view(0, 2)];
        ndarray::shallow(_kPPD) = _kFull[ndarray::view(2, _kFull.getSize<0>())];
        ndarray::Array<Pixel,3,1> finalFull(window(_xFull, _parent->_innerBBox));
        ndarray::shallow(_finalTD) = finalFull[ndarray::view(0, 2)];
        ndarray::shallow(_finalPPD) = finalFull[ndarray::view(2, finalFull.getSize<0>())];
    }

private:
    void recompute() {
        if (_valid) return;

        lsst::afw::math::FourierCutout::Ptr fourierCutout = 
            _parent->_kernelVisitor->getFourierImage();
        
        ndarray::FourierArray<Pixel, 2, 2> kernelImage(
            fourierCutout->getImageWidth(),
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
        _kTD[0] = _kTD[1] = _parent->_linearMatrixHandler->computeUnconvolvedImage() *
            kernelImage;
        ndarray::differentiate(1, _kTD[0]);
        ndarray::differentiate(0, _kTD[1]);
        _kPPD = _parent->_getMorphologyProjection()->computeProjectedParameterDerivative() *
            _parent->getPhotFactor();
        _parent->_shifter->apply(_kPPD.begin(), _kPPD.end());
        _parent->_applyKernel(_kPPD.begin(), _kPPD.end());
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

// -- multifit::FourierModelProjection::PsfMatrixHandler -----------------------

class multifit::FourierModelProjection::PsfMatrixHandler : boost::noncopyable {
public:
    void handleParameterChange() { _valid = false; }
    
    ndarray::Array<Pixel const,3,1> computePsfParameterDerivative() {
        if (!_valid) {
            lsst::afw::math::FourierCutout::Ptr fourierCutout = 
                _parent->_kernelVisitor->getFourierDerivativeImageList()[0];

            ndarray::FourierArray<Pixel, 3, 3> fourierDerivative(
                fourierCutout->getImageWidth(),
                ndarray::external(
                    fourierCutout->begin(),
                    ndarray::makeVector(
                        _parent->_kernelVisitor->getNParameters(), 
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
            _k = fourierDerivative;
            _parent->_applyKernel(
                _k.begin(), _k.end(),
                _parent->_linearMatrixHandler->computeUnconvolvedImage()
            );
            _ifft->execute();
            _valid = true;
        }
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
    //TODO: test if _kernelVisitor is null, throw exception if needed
    return _kernelVisitor->getNParameters();
}

void multifit::FourierModelProjection::_convolve(
    lsst::afw::math::Kernel::ConstPtr const & kernel
) { 
    lsst::afw::geom::PointD point = _getPsfPosition(); 
    _kernelVisitor.reset(lsst::afw::image:PointD(point.getX(), point.getY()));
    if (!_kernelVisitor->hasDerivatives()) {
        disableProducts(PSF_PARAMETER_DERIVATIVE);
    }
    if (getActiveProducts() & PSF_PARAMETER_DERIVATIVE) {
        _psfMatrixHandler.reset(new PsfMatrixHandler(this));
    } else {
        _psfMatrixHandler.reset();
    }
    _linearMatrixHandler->handleNonlinearParameterChange();
    if (_nonlinearMatrixHandler) {
        _nonlinearMatrixHandler->handleParameterChange();
    }
}

void multifit::FourierModelProjection::_computeLinearParameterDerivative(
    ndarray::Array<Pixel,2,2> const & output
) {
    _wf->compress(
        _linearMatrixHandler->computeLinearParameterDerivative(),
        output
    );
}

void multifit::FourierModelProjection::_computePsfParameterDerivative(
    ndarray::Array<Pixel,2,2> const & output
) {
    _wf->compress(
        _psfMatrixHandler->computePsfParameterDerivative(),
        output
    );
}

void multifit::FourierModelProjection::_computeTranslationDerivative(
    ndarray::Array<Pixel,2,2> const & output
) {
    _wf->compress(
        _nonlinearMatrixHandler->computeTranslationDerivative(),
        output
    );
}

void multifit::FourierModelProjection::_computeProjectedParameterDerivative(
    ndarray::Array<Pixel,2,2> const & output
) {
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
    _shifter->handleNonlinearParameterChange();
}

int multifit::FourierModelProjection::_enableProducts(int toAdd) {
    if (toAdd & PSF_PARAMETER_DERIVATIVE) {
        if (!_kernelVisitor->hasDerivative()) toAdd |= (~PSF_PARAMETER_DERIVATIVE);
    }
    if ((toAdd & NONLINEAR_PARAMETER_DERIVATIVE) && !_nonlinearMatrixHandler) {
        _nonlinearMatrixHandler.reset(new NonlinearMatrixHandler(this));
    }
    if ((toAdd & PSF_PARAMETER_DERIVATIVE) && !_psfMatrixHandler) {
        _psfMatrixHandler.reset(new PsfMatrixHandler(this));
    }
    return toAdd;
}

int multifit::FourierModelProjection::_disableProducts(int toRemove) {
    toRemove &= (~MODEL_IMAGE) & (~LINEAR_PARAMETER_DERIVATIVE);
    return toRemove;
}

multifit::FourierModelProjection::FourierModelProjection(
    ComponentModel::ConstPtr const & model,
    Kernel::ConstPtr const & kernel,
    Wcs::ConstPtr const & wcs,
    Footprint::ConstPtr const & footprint,
    double photFactor,
    int activeProducts
) : ComponentModelProjection(model,kernel,wcs,footprint,photFactor),
    _kernelVisitor(), _wf(), _outerBBox(), _innerBBox()
{
    enableProducts(activeProducts);
    _convolve(kernel);
    _setDimensions();
}

void multifit::FourierModelProjection::_setDimensions() {
    lsst::afw::geom::Extent2I dimensions = getMorphologyProjection()->getDimensions();
    lsst::afw::geom::Point2D centerOnExposure = _getPsfPosition();
    lsst::afw::geom::Point2I bboxMin = lsst::afw::geom::Point2I::makeXY(
        int(std::floor(centerOnExposure.getX() - dimensions.getX()/2)),
        int(std::floor(centerOnExposure.getY() - dimensions.getY()/2))
    );
    lsst::afw::geom::Point2I bboxMax = bboxMin + dimensions;
    _outerBBox = lsst::afw::geom::Box2I(bboxMin, bboxMax);
    // Right now, _innerBBox is defined relative to exposure    
    int padding = getMorphologyProjection()->getPadding();
    _innerBBox = _outerBBox;
    _innerBBox.expand(-padding);

    //TODO: convert footprint to use geom::Box2I
    //_innerBBox.setIntersection(getFootprint()->getBBox());
    _wf = boost::make_shared<WindowedFootprint>(*getFootprint(), _innerBBox);

    // But now, and forevermore, _innerBBox is defined relative to _outerBBox.
    _innerBBox.shift(-_outerBBox.getMin());
    _kernelVisitor->fft(_outerBBox.getWidth, _outerBBox.getHeight());    
    _linearMatrixHandler.reset(new LinearMatrixHandler(this));
    if (getActiveProducts() & NONLINEAR_PARAMETER_DERIVATIVE)
        _nonlinearMatrixHandler.reset(new NonlinearMatrixHandler(this));
    if (getActiveProducts() & PSF_PARAMETER_DERIVATIVE)
        _psfMatrixHandler.reset(new PsfMatrixHandler(this));
}

void multifit::FourierModelProjection::_applyKernel(
    ndarray::FourierArray<Pixel,3,3>::Iterator iter, 
    ndarray::FourierArray<Pixel,3,3>::Iterator const & end,
    ndarray::FourierArray<Pixel const,2,2> const & kernel
) const {
    for (; iter != end; ++iter) {
        *iter *= kernel;
    }    
}

void multifit::FourierModelProjection::_applyKernel(
    ndarray::FourierArray<Pixel,3,3>::Iterator iter, 
    ndarray::FourierArray<Pixel,3,3>::Iterator const & end
) const {
    lsst::afw::math::FourierCutout::Ptr fourierCutout = 
        _kernelVisitor->getFourierImage();
        
    ndarray::FourierArray<Pixel, 2, 2> kernelImage(
        fourierCutout->getImageWidth(),
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

    _applyKernel(iter, end, kernelImage);
}

multifit::FourierModelProjection::~FourierModelProjection() {}
