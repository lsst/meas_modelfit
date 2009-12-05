#ifndef LSST_MEAS_MULTIFIT_FOURIERMODEL_H
#define LSST_MEAS_MULTIFIT_FOURIERMODEL_H

#include <ndarray/fft.hpp>

#include <lsst/afw/math/ellipses.h>
#include <lsst/afw/math/AffineTransform.h>
#include <lsst/afw/math/ConvolutionVisitor.h>

#include <lsst/meas/multifit/Model.h>
#include <lsst/meas/multifit/ModelFactory.h>
#include <lsst/meas/multifit/ModelProjection.h>
#include <lsst/meas/multifit/WindowedFootprint.h>

#include <lsst/meas/multifit/fourier/Types.h>
#include <lsst/meas/multifit/fourier/PositionComponent.h>
#include <lsst/meas/multifit/fourier/MorphologyComponent.h>
#include <lsst/meas/multifit/fourier/Handlers.h>

namespace lsst{
namespace meas {
namespace multifit {

class FourierModelFactory;

class FourierModel : public Model {
public:
    typedef boost::shared_ptr<FourierModel> Ptr;
    typedef boost::shared_ptr<const FourierModel> CosntPtr;

    typedef fourier::PositionComponent Position;
    typedef fourier::MorphologyComponent Morphology;

    virtual Model::Ptr clone() const { return Model::Ptr(new FourierModel(*this)); }
    virtual void setNonlinearParameters(ParameterConstIterator parameterIter) {
        _position->setParameters(parameterIter);
        Model::setNonlinearParameters(parameterIter);
    }
protected:
    virtual ModelProjection::Ptr createProjection() const;
private:
    friend class FourierModelFactory;

    FourierModel(
        boost::shared_ptr<const ModelFactory> const & factory, 
        int const linearParameterSize
    ) : Model(factory, linearParameterSize) {}

    FourierModel(FourierModel const & model) : Model(model) {}


    //because position is in ra/dec, it does not need to be transformed
    //and thus can be shared between all projections
    Position::Ptr _position;

    //used for constructing projections
    Morphology::Factory::ConstPtr _morphologyFactory;
};

class FourierModelFactory : public ModelFactory {
public:
    typedef boost::shared_ptr<FourierModelFactory> Ptr;
    typedef boost::shared_ptr<FourierModelFactory const> ConstPtr;
    typedef fourier::PositionComponent Position;
    typedef fourier::MorphologyComponent Morphology;

    virtual int const getNonlinearParameterSize() const {
        return _position->getParameterSize() + 
            _morphologyFactory->getNonspatialParameterSize() +
            (_morphologyFactory->hasEllipse())? 3: 0;
    }

    virtual int const getMinLinearParameterSize() const {
        return _morphologyFactory->getMinLinearParameterSize();
    }

    virtual int const getMaxLinearParameterSize() const {
        return _morphologyFactory->getMaxLinearParameterSize();
    }

    virtual Model::Ptr makeModel(
        int linearParameterSize,
        ParameterConstIterator linearParameters,
        ParameterConstIterator nonlinearParameters
    ) const;

    static ConstPtr makeFactory(
        fourier::PositionComponent::ConstPtr const & position, 
        fourier::MorphologyComponent::Factory::ConstPtr const & morphologyFactory
    ) {
        return ConstPtr(new FourierModelFactory(position, morphologyFactory));
    }

private:

    FourierModelFactory(
        Position::ConstPtr const & position, 
        Morphology::Factory::ConstPtr const & morphologyFactory
    ) :
        _position(position),
        _morphologyFactory(morphologyFactory)
    {}

    Position::ConstPtr _position;
    Morphology::Factory::ConstPtr _morphologyFactory;
};



class FourierModelProjection :
    public ModelProjection, 
    public lsst::afw::math::Convolvable 
{
public:
    typedef boost::shared_ptr<FourierModelProjection> Ptr;
    typedef boost::shared_ptr<const FourierModelProjection> ConstPtr;

    typedef fourier::Pixel Pixel;
    typedef fourier::Complex Complex;

    typedef fourier::PositionComponent Position;
    typedef fourier::MorphologyComponent Morphology;

    static const int WCS_PARAMETER_SIZE = 6;

    virtual int const getWcsParameterSize() const;
    virtual int const getPsfParameterSize() const;

    /**
     *  Compute the model image, flattened via its footprint.
     */
    virtual void computeModelImage(ndarray::Array<Pixel,1,1> const & vector);

    /**
     *  Compute the derivative of the model image with respect to its linear 
     *  parameters, flattened via its footprint such that the first dimension 
     *  runs over parameters and the second runs over pixels.
     */
    virtual void computeLinearParameterDerivative(
        ndarray::Array<Pixel,2,2> const & matrix
    );

    /**
     *  Compute the derivative of the model image with respect to its nonlinear
     *  parameters, flattened via its footprint such that the first dimension 
     *  runs over parameters and the second runs over pixels.
     */
    virtual void computeNonlinearParameterDerivative(
        ndarray::Array<Pixel,2,2> const & matrix
    );

    /**
     *  Compute the derivative of the model image with respect to its Wcs 
     *  parameters, flattened via its footprint such that the first dimension 
     *  runs over parameters and the second runs over pixels.
     */
    virtual void computeWcsParameterDerivative(
        ndarray::Array<Pixel,2,2> const & matrix
    );

    /**
     *  Compute the derivative of the model image with respect to its PSF 
     *  parameters, flattened via its footprint such that the first dimension 
     *  runs over parameters and the second runs over pixels.
     */
    virtual void computePsfParameterDerivative(
        ndarray::Array<Pixel,2,2> const & matrix
    );

    FourierModel const & getModel() const {
        return static_cast<FourierModel const &>(*ModelProjection::getModel());
    } 

protected:

    // ------------------------- PSF/KERNEL API -------------------------- //

    /**
     *  Return the center of the object in the pixel coordinates of
     *  the Wcs it was projected onto.
     *
     *  This is the point at which the PSF model will be evaluated.
     */
    virtual lsst::afw::image::PointD getPsfPosition() const;
    virtual void setKernel(boost::shared_ptr<const Kernel> const & kernel);

    virtual void convolve(
        lsst::afw::math::FourierConvolutionVisitor const & visitor
    );

    /**
     *  Return 'true' if convolve has been called.
     */
    virtual bool isConvolved() const { return _kernelVisitor; }

    /**
     *  Return the supported convolution type flag.
     */
    virtual lsst::afw::math::ConvolutionVisitor::TypeFlag getConvolutionType() const {
        return lsst::afw::math::ConvolutionVisitor::FOURIER;
    }

    // ------------------------------------------------------------------- //

    virtual int _enableProducts(int toAdd);
    virtual int _disableProducts(int toRemove);

    virtual void _handleLinearParameterChange();
    virtual void _handleNonlinearParameterChange();

private:
    friend class FourierModel;
    friend class fourier::ModelImageHandler;
    friend class fourier::LinearMatrixHandler;
    friend class fourier::NonlinearMatrixHandler;
    friend class fourier::WcsMatrixHandler;
    friend class fourier::PsfMatrixHandler;

    typedef ndarray::FourierTransform<Pixel,2> FFT;

    FourierModelProjection(
        boost::shared_ptr<const Model> const & model,
        fourier::PositionComponent::Ptr const & position,
        fourier::MorphologyComponent::Factory::ConstPtr const & morphologyFactory
    ) : 
        ModelProjection(model),
        _position(position),
        _morphology(
            morphologyFactory->makeMorphology(getLinearParameterSize())
        )
    {
        if(_morphology->hasEllipse())
            _ellipse.reset(_morphology->getEllipse().clone());
    }
   
    void setDimensions();

    void applyShift(
        ndarray::FourierArray<Pixel,3,3>::Iterator iter, 
        ndarray::FourierArray<Pixel,3,3>::Iterator const & end
    ) const;

    void applyKernel(
        ndarray::FourierArray<Pixel,3,3>::Iterator iter, 
        ndarray::FourierArray<Pixel,3,3>::Iterator const & end,
        ndarray::FourierArray<Pixel,2,2> const & kernel
    ) const;

    void applyKernel (
        ndarray::FourierArray<Pixel,3,3>::Iterator iter, 
        ndarray::FourierArray<Pixel,3,3>::Iterator const & end
    ) const; 

    /** transform from parameter Wcs to exposure WCS */
    lsst::afw::math::AffineTransform::ConstPtr _transform;     
    /** PSF/Kernel information. */
    lsst::afw::math::FourierConvolutionVisitor::Ptr _kernelVisitor; 
    /** Determines RA/DEC center from parameters.*/
    Position::Ptr _position;     
    /** A simpler uncentered model class in the exposure Wcs. */
    Morphology::Ptr _morphology; 
    /** Size and ellipticity parameters. Empty if !_morphology->hasEllipse().*/
    lsst::afw::math::ellipses::Core::Ptr _ellipse; 
    /** maps footprint to handler output arrays */
    WindowedFootprint::ConstPtr _wf; 
    /** bounding box of padded arrays relative to exposure */
    lsst::afw::image::BBox _outer; 
    /** bounding box of unpadded arrays relative to _outer */
    lsst::afw::image::BBox _inner; 
    /** derivatives of ellipse in exposure Wcs */
    lsst::afw::math::ellipses::Core::TransformDerivative::Ptr _td; 

    boost::scoped_ptr<fourier::ModelImageHandler> _modelImageHandler;
    boost::scoped_ptr<fourier::LinearMatrixHandler> _linearMatrixHandler;
    boost::scoped_ptr<fourier::NonlinearMatrixHandler> _nonlinearMatrixHandler;
    boost::scoped_ptr<fourier::WcsMatrixHandler> _wcsMatrixHandler;
    boost::scoped_ptr<fourier::PsfMatrixHandler> _psfMatrixHandler;
};

inline lsst::afw::image::PointD FourierModelProjection::getPsfPosition() const { 
    return (*_transform)(_position->getCenter()); 
}

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_FOURIERMODEL_H
