#include <lsst/meas/multifit/FourierModel.h>
#include <lsst/afw/math/FourierCutout.h>
#include <lsst/pex/exceptions/Runtime.h>
#include <lsst/meas/multifit/Utils.h>

namespace multifit = lsst::meas::multifit;



//-------------------- FourierModel::LinearMatrixHandler --------------------//
class multifit::FourierModel::LinearMatrixHandler : boost::noncopyable {
public:
    void invalidate() { _valid = false; }
    
    int const getSize() const { return _output.size(); }

    ndarray::Array<Pixel,3,1> compute() {
        if (!_valid) 
            recompute();
        return _output; 
    }

    ndarray::FourierArray<Pixel,3,3> computeUnconvolved() {
        if (!_valid) 
            recompute();
        return _unconvolved; 
    }

    explicit LinearMatrixHandler(FourierModel * parent);

private:
    void recompute();

    bool _valid;
    FourierModel * _parent;
    FFT::Ptr _ifft;
    ndarray::Array<Pixel,3,3> _padded;
    ndarray::Array<Pixel,3,1> _output;
    ndarray::FourierArray<Pixel,3,3> _unconvolved;
    ndarray::FourierArray<Pixel,3,3> _workspace;
};



multifit::FourierModel::LinearMatrixHandler::LinearMatrixHandler(
    multifit::FourierModel * parent
) : 
    _valid(false), _parent(parent) 
{
    ndarray::Vector<int,3> shape = ndarray::makeVector(
        _parent->_morphology->getLinearParameterSize(),
        _parent->_outer.getHeight(),
        _parent->_outer.getWidth()
    );
    _ifft = FFT::planMultiplexInverse(shape,_workspace,_padded);
    ndarray::shallow(_unconvolved) = FFT::initializeK(shape);
    ndarray::shallow(_output) = window(_padded,_parent->_inner);
}

void multifit::FourierModel::LinearMatrixHandler::recompute() {
    _workspace = _parent->_morphology->computeLinearParameterDerivative();
    _workspace *= _parent->getPhotFactor();
    _parent->applyShift(_workspace.begin(),_workspace.end());
    _unconvolved = _workspace;
    _parent->applyKernel(_workspace.begin(),_workspace.end());
    _ifft->execute();
    _valid = true;
}

//--------------------- FourierModel::ModelImageHandler ---------------------//
class multifit::FourierModel::ModelImageHandler : boost::noncopyable {
public:
    void invalidate() { _valid = false; }

    ndarray::FourierArray<Pixel,2,2> computeUnconvolved() {
        if (!_valid) 
            recompute();
        return _unconvolved;
    }
    ndarray::FourierArray<Pixel,3,3> computeConvolvedDXY() {
        if (!_valid) 
            recompute();
        return _convolvedDXY; 
    }

    explicit ModelImageHandler(FourierModel * parent);

private:
    void recompute();

    bool _valid;
    FourierModel * _parent;
    ndarray::FourierArray<Pixel,2,2> _unconvolved;
    ndarray::FourierArray<Pixel,3,3> _convolvedDXY;
};

multifit::FourierModel::ModelImageHandler::ModelImageHandler(
    multifit::FourierModel * parent
) :
    _valid(false), _parent(parent)
{
    ndarray::Vector<int,2> shape = ndarray::makeVector(
        _parent->_outer.getHeight(),
        _parent->_outer.getWidth()
    );
    ndarray::shallow(_unconvolved) = FFT::initializeK(shape);
    ndarray::shallow(_convolvedDXY) = FFT::initializeK(
        ndarray::concatenate(2,shape)
    );
}

void multifit::FourierModel::ModelImageHandler::recompute() {
    lsst::afw::math::FourierCutout::Ptr cutout(
        _parent->_kernelVisitor->getFourierImage()
    );
    boost::shared_ptr<Complex> owner(cutout->getOwner());
    int width = cutout->getFourierWidth();
    int height = cutout->getFourierHeight();
    ndarray::FourierArray<Pixel, 2, 2> parentKernelImage(
        cutout->getImageWidth(),
        ndarray::external(
            owner.get(),
            ndarray::makeVector(height, width),
            ndarray::makeVector(width, 1),
            owner
        )
    );

    ndarray::FourierArray<Pixel,3,3> linearMatrix = 
            _parent->_linearMatrixHandler->computeUnconvolved();
    getVectorView(_unconvolved.getBase()) = 
            getMatrixView(linearMatrix.getBase()) 
            * _parent->getLinearParameters();
    _convolvedDXY[0] = _convolvedDXY[1] = _unconvolved * parentKernelImage;
    ndarray::differentiate(0, _convolvedDXY[1]);
    ndarray::differentiate(1, _convolvedDXY[0]);
    _valid = true;
}

// ------------------ FourierModel::NonlinearMatrixHandler ------------------ //

class multifit::FourierModel::NonlinearMatrixHandler : boost::noncopyable {
public:
    void invalidate() { _valid = false; }

    int const getSize() const { return _output.getSize<0>(); }

    ndarray::Array<Pixel,3,1> compute() {
        if (!_valid) recompute();
        return _output;
    }

    explicit NonlinearMatrixHandler(FourierModel * parent);

private:
    void recompute();

    bool _valid;
    FourierModel * _parent;
    FFT::Ptr _ifft;
    ndarray::Array<Pixel,3,3> _padded;
    ndarray::Array<Pixel,3,1> _output;
    ndarray::FourierArray<Pixel,3,3> _workspace;
    ndarray::FourierArray<Pixel,3,3> _dPosition;
    ndarray::FourierArray<Pixel,3,3> _dEllipse;
    ndarray::FourierArray<Pixel,3,3> _dNonspatial;
};

multifit::FourierModel::NonlinearMatrixHandler::NonlinearMatrixHandler(
        multifit::FourierModel * parent
) :
    _valid(false), _parent(parent)
{
    int last = _parent->_position->getParameterSize();
    int size = last + _parent->_morphology->getNonspatialParameterSize();
    if (_parent->_ellipse) 
        size += 3;
    ndarray::Vector<int,3> shape = ndarray::makeVector(
        size,
        _parent->_outer.getHeight(),
        _parent->_outer.getWidth()
    );

    _ifft = FFT::planMultiplexInverse(shape,_workspace,_padded);
    ndarray::shallow(_output) = window(_padded,_parent->_inner);

    ndarray::shallow(_dPosition) = _workspace[ndarray::view(0,last)];
    if (_parent->_ellipse) {
        ndarray::shallow(_dEllipse) = _workspace[ndarray::view(last,last+3)];
        last += 3;
    }
    if (_parent->_morphology->getNonspatialParameterSize() > 0) {
        ndarray::shallow(_dNonspatial) = _workspace[ndarray::view(last,size)];
    }
}

void multifit::FourierModel::NonlinearMatrixHandler::recompute() {
    ndarray::FourierArray<Pixel,3,3> convolvedDXY( 
            _parent->_modelImageHandler->computeConvolvedDXY()
    );
    getMatrixView(_dPosition.getBase()) = 
            getMatrixView(convolvedDXY.getBase())
            * _parent->_transform->matrix().linear()
            * _parent->_position->getCenterDerivative();
    if (_parent->_ellipse) {        
        getMatrixView(_dEllipse.getBase()) = 
            getMatrixView(
                _parent->_morphology->computeEllipseDerivative().getBase()
            ) * _parent->_td->dInput();    
    }
    if (_parent->_morphology->getNonspatialParameterSize() > 0) {
        _dNonspatial = 
                _parent->_morphology->computeNonspatialParameterDerivative();
    }
    _parent->applyShift(
            _workspace.begin()+_parent->_position->getParameterSize(),
            _workspace.end()
    );
    _parent->applyKernel(
        _workspace.begin()+_parent->_position->getParameterSize(),
        _workspace.end()
    );
    _ifft->execute();
    _valid = true;
}

// --------------------- FourierModel::WcsMatrixHandler --------------------- //

class multifit::FourierModel::WcsMatrixHandler : boost::noncopyable {
public:
    void invalidate() { _valid = false; }

    int const getSize() const { return _output.getSize<0>(); }

    ndarray::Array<Pixel,3,1> compute() {
        if (!_valid) recompute();
        return _output;
    }

    explicit WcsMatrixHandler(FourierModel * parent);

private:
    void recompute();

    bool _valid;
    FourierModel * _parent;
    FFT::Ptr _ifft;
    ndarray::Array<Pixel,3,3> _padded;
    ndarray::Array<Pixel,3,1> _output;
    ndarray::FourierArray<Pixel,3,3> _workspace;
};

multifit::FourierModel::WcsMatrixHandler::WcsMatrixHandler(
        FourierModel * parent
) :
    _valid(false), _parent(parent)
{
    _ifft = FFT::planMultiplexInverse(
        ndarray::makeVector(
            WCS_PARAMETER_SIZE,
            _parent->_outer.getHeight(),
            _parent->_outer.getWidth()
        ),
        _workspace,
        _padded
    );
    ndarray::shallow(_output) = window(_padded,_parent->_inner);
}

void multifit::FourierModel::WcsMatrixHandler::recompute() {
    Eigen::Map<Eigen::Matrix<Complex,Eigen::Dynamic,Eigen::Dynamic> > workspaceView = getMatrixView(
        _workspace.getBase()
    );
    workspaceView = 
        getMatrixView(_parent->_morphology->computeEllipseDerivative().getBase()) 
        * _parent->_td->dTransform();
    _parent->applyShift(_workspace.begin(), _workspace.end());
    _parent->applyKernel(_workspace.begin(),_workspace.end());
    ndarray::FourierArray<Pixel,3,3> convolvedDXY = 
        _parent->_modelImageHandler->computeConvolvedDXY();
    workspaceView += getMatrixView(convolvedDXY.getBase()) 
        * _parent->_transform->d(_parent->_position->getCenter());
    
    _ifft->execute();
    _valid = true;
}

// --------------------- FourierModel::PsfMatrixHandler --------------------- //

class multifit::FourierModel::PsfMatrixHandler : boost::noncopyable {
public:
    void invalidate() { _valid = false; }
    
    int const getSize() const { return _output.getSize<0>(); }

    ndarray::Array<Pixel,3,1> compute() {
        if (!_valid) recompute();
        return _output;
    }

    explicit PsfMatrixHandler(FourierModel * parent);

private:
    void recompute();

    bool _valid;
    FourierModel * _parent;
    FFT::Ptr _ifft;
    ndarray::Array<Pixel,3,3> _padded;
    ndarray::Array<Pixel,3,1> _output;
    ndarray::FourierArray<Pixel,3,3> _workspace;
};

multifit::FourierModel::PsfMatrixHandler::PsfMatrixHandler(
        FourierModel * parent
) :
    _valid(false), _parent(parent)
{
    _ifft = FFT::planMultiplexInverse(
        ndarray::makeVector(
            _parent->_kernelVisitor->getNParameters(),
            _parent->_outer.getHeight(),
            _parent->_outer.getWidth()
        ),
        _workspace,
        _padded
    );
    ndarray::shallow(_output) = window(_padded,_parent->_inner);
}

void multifit::FourierModel::PsfMatrixHandler::recompute() {
    ndarray::FourierArray<Pixel,2,2> unconvolved( 
            _parent->_modelImageHandler->computeUnconvolved()
    );
    std::vector<lsst::afw::math::FourierCutout::Ptr> cutoutList = 
        _parent->_kernelVisitor->getFourierDerivativeImageList();
    lsst::afw::math::FourierCutout::Ptr first = cutoutList[0];
    boost::shared_ptr<Complex> owner(first->getOwner());
    int width = first->getFourierWidth();
    int height = first->getFourierHeight();
    int nDerivatives = static_cast<int>(cutoutList.size());
    _workspace = ndarray::FourierArray<Pixel, 3, 3>(
        first->getImageWidth(),
        ndarray::external(
            owner.get(),
            ndarray::makeVector(nDerivatives, height, width),
            ndarray::makeVector(height*width, width, 1),
            owner
        )
    );
       
    _parent->applyKernel(_workspace.begin(),_workspace.end(),unconvolved);
    _ifft->execute();
    _valid = true;
}

// ----------------------------- FourierModel ------------------------------- //

void multifit::FourierModel::reproject(
    Kernel const & kernel,
    WcsConstPtr const & wcs,
    FootprintConstPtr const & footprint,
    double photFactor
) {
    setProjectionVariables(wcs,footprint,photFactor);
    convolve(*kernel.computeFourierConvolutionVisitor(getPsfPosition()));
}

void multifit::FourierModel::_handleLinearParameterChange() {
    _morphology->setLinearParameters(getLinearParameters().data());
    if (_morphology->getDimensions() != _outer.getDimensions()) {
        setDimensions();
    } else {
        _modelImageHandler->invalidate();
        if (_nonlinearMatrixHandler) _nonlinearMatrixHandler->invalidate();
        if (_wcsMatrixHandler) _wcsMatrixHandler->invalidate();
        if (_psfMatrixHandler) _psfMatrixHandler->invalidate();
    }
}

void multifit::FourierModel::_handleNonlinearParameterChange() {
    ParameterConstIterator parameters = getNonlinearParameters().data();
    _position->setParameters(parameters);
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

void multifit::FourierModel::computeModelImage(
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
    v = m * getLinearParameters();
}

void multifit::FourierModel::computeLinearParameterDerivative(
        ndarray::Array<Pixel,2,2> const & output
) {
    assert(output.getSize<0>() != getLinearParameterSize());
    assert(output.getSize<1>() == _wf->getNpix());
    _wf->compress(_linearMatrixHandler->compute(),output);
}

void multifit::FourierModel::computeNonlinearParameterDerivative(
        ndarray::Array<Pixel,2,2> const & output
) {
    assert(output.getSize<0>() != getNonlinearParameterSize());
    assert(output.getSize<1>() == _wf->getNpix());
    if (!(getActiveProducts() & NONLINEAR_PARAMETER_DERIVATIVE)) 
        return;
    _wf->compress(_nonlinearMatrixHandler->compute(),output);
}

void multifit::FourierModel::computeWcsParameterDerivative(
        ndarray::Array<Pixel,2,2> const & output
) {
    assert(output.getSize<0>() != getWcsParameterSize());
    assert(output.getSize<1>() == _wf->getNpix());
    if (!(getActiveProducts() & WCS_PARAMETER_DERIVATIVE)) 
        return;
    _wf->compress(_wcsMatrixHandler->compute(),output);
}

void multifit::FourierModel::computePsfParameterDerivative(
    ndarray::Array<Pixel,2,2> const & output
) {
    assert(output.getSize<0>() != getPsfParameterSize());
    assert(output.getSize<1>() == _wf->getNpix());
    if (!(getActiveProducts() & PSF_PARAMETER_DERIVATIVE)) 
        return;
    _wf->compress(_psfMatrixHandler->compute(),output);
}

int const multifit::FourierModel::getWcsParameterSize() const {
    if(getActiveProducts() & WCS_PARAMETER_DERIVATIVE != 0) 
        return _wcsMatrixHandler->getSize();
    else return 0;
}

int const multifit::FourierModel::getPsfParameterSize() const {
    if(getActiveProducts() & PSF_PARAMETER_DERIVATIVE != 0) 
        return _psfMatrixHandler->getSize();
    else return 0;
}

void multifit::FourierModel::convolve(
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
        _psfMatrixHandler.reset(new PsfMatrixHandler(this));
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

int multifit::FourierModel::_enableProducts(int toAdd) {
    if (toAdd & PSF_PARAMETER_DERIVATIVE) {
        if (!_kernelVisitor->hasDerivatives()) 
            toAdd |= (~PSF_PARAMETER_DERIVATIVE);
    }
    if ((toAdd & NONLINEAR_PARAMETER_DERIVATIVE) && !_nonlinearMatrixHandler) {
        _nonlinearMatrixHandler.reset(new NonlinearMatrixHandler(this));
    }
    if ((toAdd & WCS_PARAMETER_DERIVATIVE) && !_wcsMatrixHandler) {
        _wcsMatrixHandler.reset(new WcsMatrixHandler(this));
    }
    if ((toAdd & PSF_PARAMETER_DERIVATIVE) && !_psfMatrixHandler) {
        _psfMatrixHandler.reset(new PsfMatrixHandler(this));
    }
    return toAdd;
}

int multifit::FourierModel::_disableProducts(int toRemove) {
    toRemove &= (~MODEL_IMAGE) & (~LINEAR_PARAMETER_DERIVATIVE);
    return toRemove;
}

void multifit::FourierModel::setDimensions() {
    std::pair<int, int> dimensions = _morphology->getDimensions();
    //TODO: add method to compute linear approximation of a WCS
    //at a given position:

    //replace following line with a call to that method
    //
    _transform.reset(new lsst::afw::math::AffineTransform(_position->getCenter()));
                     
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
   
    _linearMatrixHandler.reset(new LinearMatrixHandler(this));
    _modelImageHandler.reset(new ModelImageHandler(this));
    if (getActiveProducts() & NONLINEAR_PARAMETER_DERIVATIVE)
        _nonlinearMatrixHandler.reset(new NonlinearMatrixHandler(this));
    if (getActiveProducts() & WCS_PARAMETER_DERIVATIVE)
        _wcsMatrixHandler.reset(new WcsMatrixHandler(this));
    if (getActiveProducts() & PSF_PARAMETER_DERIVATIVE)
        _psfMatrixHandler.reset(new PsfMatrixHandler(this));
}

void multifit::FourierModel::applyShift(
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

void multifit::FourierModel::applyKernel(
    ndarray::FourierArray<Pixel,3,3>::Iterator iter, 
    ndarray::FourierArray<Pixel,3,3>::Iterator const & end,
    ndarray::FourierArray<Pixel,2,2> const & kernel
) const {
    for (; iter != end; ++iter) {
        *iter *= kernel;
    }    
}

void multifit::FourierModel::applyKernel(
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


multifit::FourierModel::FourierModel(
    FourierModel const & other
) : 
    Model(other),
    _transform(other._transform),
    _position(other._position->clone()),
    _morphology(other._morphology->clone()),
    _ellipse(other._ellipse->clone()),
    _wf(other._wf),
    _outer(other._outer),
    _inner(other._inner),
    _td(other._td)
{
    convolve(*other._kernelVisitor);
}

multifit::Model::Definition multifit::FourierModel::makeDefinition(
    Position::ConstPtr const & position,
    Morphology::Factory::ConstPtr const & morphologyFactory,
    ParameterConstIterator linearBegin,
    ParameterConstIterator const linearEnd,
    ParameterConstIterator nonlinearBegin
) {
    Model::Factory::ConstPtr factory(new FourierModel::Factory(position,morphologyFactory));
    return Model::Definition(factory,linearBegin,linearEnd,nonlinearBegin);
}

multifit::FourierModel::FourierModel(
    Position::Ptr const & position, 
    Morphology::Factory::ConstPtr const & morphologyFactory,
    ParameterConstIterator linearBegin,
    ParameterConstIterator const linearEnd,
    ParameterConstIterator nonlinearBegin,
    Kernel const & kernel,
    WcsConstPtr const & wcs,
    FootprintConstPtr const & footprint,
    double photFactor
) : 
    Model(
        makeDefinition(
            position,morphologyFactory,linearBegin,linearEnd,nonlinearBegin
        ),
        MODEL_IMAGE | LINEAR_PARAMETER_DERIVATIVE,
        wcs, footprint, photFactor
    ),
    _position(position), 
    _morphology(morphologyFactory->create(linearEnd-linearBegin))
{
    convolve(*kernel.computeFourierConvolutionVisitor(getPsfPosition()));
    if (_morphology->hasEllipse()) 
        _ellipse.reset(_morphology->getEllipse().clone());
    _morphology->setLinearParameters(linearBegin);
    _handleNonlinearParameterChange();
}
    
// ------------------------- FourierModel::Factory -------------------------- //


multifit::FourierModel * multifit::FourierModel::Factory::project(
    ParameterConstIterator linearParameterBegin,
    ParameterConstIterator const linearParameterEnd,
    ParameterConstIterator nonlinearParameterBegin,
    Kernel const & kernel,
    WcsConstPtr const & wcs,
    FootprintConstPtr const & footprint,
    double photFactor
) const {
    Position::Ptr position(_position->clone());
    return new FourierModel(
        position,_morphologyFactory,
        linearParameterBegin,linearParameterEnd,
        nonlinearParameterBegin,
        kernel, wcs, footprint, photFactor
    );
}

int const multifit::FourierModel::Factory::getNonlinearParameterSize() const {
    int size = _position->getParameterSize()
            + _morphologyFactory->getNonspatialParameterSize();

    if (_morphologyFactory->hasEllipse()) 
        size += 3;
    return size;
}
