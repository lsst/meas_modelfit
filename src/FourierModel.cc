#include <lsst/meas/multifit/FourierModel.h>
#include <lsst/afw/math/FourierCutout.h>
#include <lsst/pex/exceptions/Runtime.h>
#include <lsst/meas/multifit/Utils.h>

namespace multifit = lsst::meas::multifit;

void multifit::FourierModelProjection::_handleLinearParameterChange() {
    _morphology->setLinearParameters(getModel().getLinearParameters().data());
    if (_morphology->getDimensions() != _outer.getDimensions()) {
        setDimensions();
    } else {
        _modelImageHandler->invalidate();
        if (_nonlinearMatrixHandler) _nonlinearMatrixHandler->invalidate();
        if (_wcsMatrixHandler) _wcsMatrixHandler->invalidate();
        if (_psfMatrixHandler) _psfMatrixHandler->invalidate();
    }
}

void multifit::FourierModelProjection::_handleNonlinearParameterChange() {
    ParameterConstIterator parameters = 
        getModel().getNonlinearParameters().data();
    parameters += _position->getParameterSize();
    if (_ellipse) {
        (*_ellipse)[0] = parameters[0];
        (*_ellipse)[1] = parameters[1];
        (*_ellipse)[2] = parameters[2]; 
        _ellipse->transform(*_transform);
        _morphology->setEllipse(*_ellipse);
        (*_ellipse)[0] = parameters[0];
        (*_ellipse)[1] = parameters[1];
        (*_ellipse)[2] = parameters[2];
        parameters += 3;
        *_td = lsst::afw::math::ellipses::Core::TransformDerivative(
            *_ellipse,*_transform
        );
    }
    _morphology->setNonspatialParameters(parameters);
    if (_morphology->getDimensions() != _outer.getDimensions()) {
        setDimensions();
    } else {
        _modelImageHandler->invalidate();
        _linearMatrixHandler->invalidate();
        if (_nonlinearMatrixHandler) _nonlinearMatrixHandler->invalidate();
        if (_wcsMatrixHandler) _wcsMatrixHandler->invalidate();
        if (_psfMatrixHandler) _psfMatrixHandler->invalidate();
    }
}

void multifit::FourierModelProjection::computeModelImage(
        ndarray::Array<Pixel,1,1> const & output
) {
    assert(output.getSize<0>() == _wf->getNpix());
    ndarray::Array<Pixel,2,2> compressedLinearMatrix = ndarray::allocate(
        ndarray::makeVector(getLinearParameterSize(),output.getSize<0>())
    );
    _wf->compress(_linearMatrixHandler->compute(), compressedLinearMatrix);
    Eigen::Map<Eigen::MatrixXd> m(
        compressedLinearMatrix.getData(),
        getLinearParameterSize(),
        output.getSize<0>()
    );
    Eigen::Map<Eigen::VectorXd> v(
        output.getData(),
        output.getSize<0>()
    );
    v = m * getModel().getLinearParameters();
}

void multifit::FourierModelProjection::computeLinearParameterDerivative(
        ndarray::Array<Pixel,2,2> const & output
) {
    assert(output.getSize<0>() != getLinearParameterSize());
    assert(output.getSize<1>() == _wf->getNpix());
    _wf->compress(_linearMatrixHandler->compute(),output);
}

void multifit::FourierModelProjection::computeNonlinearParameterDerivative(
        ndarray::Array<Pixel,2,2> const & output
) {
    assert(output.getSize<0>() != getNonlinearParameterSize());
    assert(output.getSize<1>() == _wf->getNpix());
    if (!(getActiveProducts() & NONLINEAR_PARAMETER_DERIVATIVE)) 
        return;
    _wf->compress(_nonlinearMatrixHandler->compute(),output);
}

void multifit::FourierModelProjection::computeWcsParameterDerivative(
        ndarray::Array<Pixel,2,2> const & output
) {
    assert(output.getSize<0>() != getWcsParameterSize());
    assert(output.getSize<1>() == _wf->getNpix());
    if (!(getActiveProducts() & WCS_PARAMETER_DERIVATIVE)) 
        return;
    _wf->compress(_wcsMatrixHandler->compute(),output);
}

void multifit::FourierModelProjection::computePsfParameterDerivative(
    ndarray::Array<Pixel,2,2> const & output
) {
    assert(output.getSize<0>() != getPsfParameterSize());
    assert(output.getSize<1>() == _wf->getNpix());
    if (!(getActiveProducts() & PSF_PARAMETER_DERIVATIVE)) 
        return;
    _wf->compress(_psfMatrixHandler->compute(),output);
}

int const multifit::FourierModelProjection::getWcsParameterSize() const {
    if(getActiveProducts() & WCS_PARAMETER_DERIVATIVE != 0) 
        return _wcsMatrixHandler->getSize();
    else return 0;
}

int const multifit::FourierModelProjection::getPsfParameterSize() const {
    if(getActiveProducts() & PSF_PARAMETER_DERIVATIVE != 0) 
        return _psfMatrixHandler->getSize();
    else return 0;
}

void multifit::FourierModelProjection::setKernel(
    boost::shared_ptr<const lsst::afw::math::Kernel> const & kernel
) {
    convolve(*kernel->computeFourierConvolutionVisitor(getPsfPosition()));
}

void multifit::FourierModelProjection::convolve(
    lsst::afw::math::FourierConvolutionVisitor const & visitor
) {
    _kernelVisitor.reset(
        new lsst::afw::math::FourierConvolutionVisitor(visitor)
    );

    int kernelWidth = lsst::afw::math::FourierCutout::computeFourierWidth(
        visitor.getWidth()
    );
    int kernelSize = visitor.getHeight() * kernelWidth;
    if (kernelSize > _morphology->getMaxKernelSize())
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Kernel exceeds model padding."
        );
    if (!visitor.hasDerivatives()) 
        disableProducts(PSF_PARAMETER_DERIVATIVE);

    if (getActiveProducts() & PSF_PARAMETER_DERIVATIVE) {
        _psfMatrixHandler.reset(new fourier::PsfMatrixHandler(this));
    } else {
        _psfMatrixHandler.reset();
    }
    _modelImageHandler->invalidate();
    _linearMatrixHandler->invalidate();
    if (_nonlinearMatrixHandler) 
        _nonlinearMatrixHandler->invalidate();
    if (_wcsMatrixHandler) 
        _wcsMatrixHandler->invalidate();
}

int multifit::FourierModelProjection::_enableProducts(int toAdd) {
    if (toAdd & PSF_PARAMETER_DERIVATIVE) {
        if (!_kernelVisitor->hasDerivatives()) 
            toAdd |= (~PSF_PARAMETER_DERIVATIVE);
    }
    if ((toAdd & NONLINEAR_PARAMETER_DERIVATIVE) && !_nonlinearMatrixHandler) {
        _nonlinearMatrixHandler.reset(
            new fourier::NonlinearMatrixHandler(this)
        );
    }
    if ((toAdd & WCS_PARAMETER_DERIVATIVE) && !_wcsMatrixHandler) {
        _wcsMatrixHandler.reset(new fourier::WcsMatrixHandler(this));
    }
    if ((toAdd & PSF_PARAMETER_DERIVATIVE) && !_psfMatrixHandler) {
        _psfMatrixHandler.reset(new fourier::PsfMatrixHandler(this));
    }
    return toAdd;
}

int multifit::FourierModelProjection::_disableProducts(int toRemove) {
    toRemove &= (~MODEL_IMAGE) & (~LINEAR_PARAMETER_DERIVATIVE);
    return toRemove;
}

void multifit::FourierModelProjection::setDimensions() {
    std::pair<int, int> dimensions = _morphology->getDimensions();
    //TODO: add method to compute linear approximation of a WCS
    //at a given position, then replace following line with a call to that 
    //function    
    _transform.reset(
        new lsst::afw::math::AffineTransform(_position->getCenter())
    );
                     
    lsst::afw::image::PointD centerOnExposure = (*_transform)(
        _position->getCenter()
    );
    lsst::afw::image::PointI offset(
        int(std::floor(centerOnExposure.getX() - dimensions.first/2)),
        int(std::floor(centerOnExposure.getY() - dimensions.second/2))
    );
    _outer = lsst::afw::image::BBox(
        offset, 
        dimensions.first, 
        dimensions.second
    );
    int padding = _morphology->getPadding();
    offset.shift(padding, padding);    
    _inner = lsst::afw::image::BBox(
        offset,
        dimensions.first - 2*padding,
        dimensions.second - 2*padding
    );
    _inner.clip(getFootprint()->getBBox());
    //_wf = window(*getFootprint(),_inner);
    
    _inner.shift(-_outer.getX0(), -_outer.getY0());
    _kernelVisitor->fft(_outer.getWidth(), _outer.getHeight());
   
    _linearMatrixHandler.reset(new fourier::LinearMatrixHandler(this));
    _modelImageHandler.reset(new fourier::ModelImageHandler(this));
    if (getActiveProducts() & NONLINEAR_PARAMETER_DERIVATIVE) {
        _nonlinearMatrixHandler.reset(
            new fourier::NonlinearMatrixHandler(this)
        );
    }
    if (getActiveProducts() & WCS_PARAMETER_DERIVATIVE)
        _wcsMatrixHandler.reset(new fourier::WcsMatrixHandler(this));
    if (getActiveProducts() & PSF_PARAMETER_DERIVATIVE)
        _psfMatrixHandler.reset(new fourier::PsfMatrixHandler(this));
}

void multifit::FourierModelProjection::applyShift(
    ndarray::FourierArray<Pixel,3,3>::Iterator iter, 
    ndarray::FourierArray<Pixel,3,3>::Iterator const & end
) const {
    lsst::afw::image::PointD center = (*_transform)(_position->getCenter());
    double dx = center.getX() - _outer.getX0();
    double dy = center.getY() - _outer.getY0();    
    for (; iter != end; ++iter) {
        ndarray::shift(
            ndarray::makeVector(dy, dx),
            *iter
        );
    }
}

void multifit::FourierModelProjection::applyKernel(
    ndarray::FourierArray<Pixel,3,3>::Iterator iter, 
    ndarray::FourierArray<Pixel,3,3>::Iterator const & end,
    ndarray::FourierArray<Pixel,2,2> const & kernel
) const {
    for (; iter != end; ++iter) {
        *iter *= kernel;
    }    
}

void multifit::FourierModelProjection::applyKernel(
    ndarray::FourierArray<Pixel,3,3>::Iterator iter, 
    ndarray::FourierArray<Pixel,3,3>::Iterator const & end
) const {
    lsst::afw::math::FourierCutout::Ptr cutout(
        _kernelVisitor->getFourierImage()
    );
    boost::shared_ptr<Complex> owner(cutout->getOwner());
    int width =cutout->getFourierWidth();
    int height = cutout->getFourierHeight();
    ndarray::FourierArray<Pixel, 2, 2> kernelImage(
        cutout->getImageWidth(),
        ndarray::external(
            owner.get(),
            ndarray::makeVector(height, width),
            ndarray::makeVector(width, 1),
            owner
        )
    );
    applyKernel(iter,end,kernelImage);
}

multifit::ModelProjection::Ptr multifit::FourierModel::createProjection() const {
    Model::ConstPtr const model(shared_from_this()); 
    return boost::shared_ptr<ModelProjection>(
        new FourierModelProjection(model, _position, _morphologyFactory)
    );
}

multifit::Model::Ptr multifit::FourierModelFactory::makeModel(
    int const linearParameterSize,
    ParameterConstIterator linearParameters,
    ParameterConstIterator nonlinearParameters
) const {
    if(linearParameterSize < getMinLinearParameterSize() || 
       linearParameterSize < getMaxLinearParameterSize()) {
        throw LSST_EXCEPT(lsst::pex::exceptions::LengthErrorException,
            (   
                boost::format("linearParameterSize must be in range [%1%, %2%]. Was %3%")
                % getMinLinearParameterSize() 
                % getMaxLinearParameterSize() 
                % linearParameterSize
             ).str()
        );
    }
    ModelFactory::ConstPtr const factory(shared_from_this());
    FourierModel::Ptr model(
        new FourierModel(
            factory ,
            linearParameterSize
        )
    );
    model->_position.reset(_position->clone());
    model->_morphologyFactory = _morphologyFactory;

    model->setLinearParameters(linearParameters);
    model->setNonlinearParameters(nonlinearParameters);

    return model;
}


