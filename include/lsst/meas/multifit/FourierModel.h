#ifndef LSST_MEAS_MULTIFIT_FOURIERMODEL_H
#define LSST_MEAS_MULTIFIT_FOURIERMODEL_H

#include <ndarray/fft.hpp>

#include <lsst/afw/math/ellipses.h>
#include <lsst/meas/multifit/Model.h>
#include <lsst/meas/multifit/WindowedFootprint.h>
#include <lsst/afw/math/AffineTransform.h>
#include <lsst/afw/math/ConvolutionVisitor.h>

namespace lsst{
namespace meas {
namespace multifit {

class FourierModel : public Model, public lsst::afw::math::Convolvable {
public:
    typedef boost::shared_ptr<FourierModel> Ptr;
    typedef boost::shared_ptr<const FourierModel> ConstPtr;

    typedef Model::Pixel Pixel;
    typedef std::complex<Pixel> Complex;
    typedef ndarray::FourierTransform<Pixel,2> FFT;

    class Factory;
    class Position;
    class Morphology;

    static const int WCS_PARAMETER_SIZE = 6;

    virtual FourierModel * clone() const { return new FourierModel(*this); }

    virtual void reproject(
        Kernel const & kernel,
        WcsConstPtr const &wcs,
        FootprintConstPtr const &footprint,
        double photFactor
    );

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

    class Morphology {
    public:
        typedef boost::shared_ptr<Morphology> Ptr;
        typedef boost::shared_ptr<const Morphology> ConstPtr;

        class Factory {
        public:
            typedef boost::shared_ptr<Factory> Ptr;
            typedef boost::shared_ptr<const Factory> ConstPtr;

            virtual Morphology * create(int linearParameterSize) const = 0;

            virtual bool hasEllipse() const = 0;

            virtual int const getMinLinearParameterSize() const = 0;
            virtual int const getMaxLinearParameterSize() const = 0;
            virtual int const getNonspatialParameterSize() const = 0;

            virtual ~Factory() {}
        };

        virtual Morphology * clone() const = 0;

        virtual void getLinearParameters(
            ParameterIterator parameters
        ) const = 0;
        virtual void setLinearParameters(ParameterConstIterator parameters) = 0;

        virtual void getNonspatialParameters(
            ParameterIterator parameters
        ) const = 0;
        virtual void setNonspatialParameters(
            ParameterConstIterator parameters
        ) = 0;

        virtual lsst::afw::math::ellipses::Core const & getEllipse() const = 0;
        virtual void setEllipse(
            lsst::afw::math::ellipses::Core const & ellipse
        ) = 0;

        virtual bool hasEllipse() const = 0;
        virtual int const getLinearParameterSize() const = 0;
        virtual int const getNonspatialParameterSize() const = 0;

        //get (width, height) as a std::pair
        virtual std::pair<int, int> getDimensions() const = 0;

        virtual int const getMaxKernelSize() const = 0;
        virtual int const getPadding() const { 
            return getMaxKernelSize() / 2 + 1; 
        }

        /**
         *  Compute the derivative of the Fourier-space model with respect to 
         *  the linear parameters.
         */
        virtual ndarray::FourierArray<Pixel,3,3> computeLinearParameterDerivative() = 0;

        /**
         *  Compute the derivative of the Fourier-space model with respect to 
         *  the non-elllipse nonlinear parameters.
         */
        virtual ndarray::FourierArray<Pixel,3,3> computeNonspatialParameterDerivative() = 0;

        /**
         *  Compute the derivative of the Fourier-space model with respect to 
         *  the ellipse parameters.
         */
        virtual ndarray::FourierArray<Pixel,3,3> computeEllipseDerivative() = 0;

        virtual ~Morphology() {}

    private:
        void operator=(Morphology const & other) {}
    };

    class Position {
    public:
        typedef boost::shared_ptr<Position> Ptr;
        typedef boost::shared_ptr<const Position> ConstPtr;
        typedef Eigen::Matrix<double, 2, Eigen::Dynamic> CenterDerivative;
        virtual ~Position() {}

        virtual Position * clone() const = 0;

        virtual void getParameters(ParameterIterator parameters) const = 0;
        virtual void setParameters(ParameterConstIterator parameters) = 0;

        virtual int const getParameterSize() const = 0;

        /**
         *  Return the center of the object in (ra,dec).
         */
        virtual lsst::afw::image::PointD const & getCenter() const = 0;

        /**
         *  Return the derivative of getCenter() with respect to the parameters.
         */
        virtual CenterDerivative const & getCenterDerivative() const = 0;

    private:
        void operator=(Position const & other) {}
    };

protected:

    // ------------------------- PSF/KERNEL API -------------------------- //

    /**
     *  Return the center of the object in the pixel coordinates of
     *  the Wcs it was projected onto.
     *
     *  This is the point at which the PSF model will be evaluated.
     */
    virtual lsst::afw::image::PointD getPsfPosition() const;

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

    FourierModel(FourierModel const & other);
    
    FourierModel(
        Position::Ptr const & position, 
        Morphology::Factory::ConstPtr const & morphologyFactory,
        ParameterConstIterator linearBegin,
        ParameterConstIterator const linearEnd,
        ParameterConstIterator nonlinearBegin,
        Kernel const & kernel,
        WcsConstPtr const &wcs,
        FootprintConstPtr const &footprint,
        double photFactor
    );  
    
    static Definition makeDefinition(
        Position::ConstPtr const & position,
        Morphology::Factory::ConstPtr const & morphologyFactory,
        ParameterConstIterator linearBegin,
        ParameterConstIterator const linearEnd,
        ParameterConstIterator nonlinearBegin
    );

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
    boost::shared_ptr<Position> _position;     
    /** A simpler uncentered model class in the exposure Wcs. */
    boost::shared_ptr<Morphology> _morphology; 
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

    class ModelImageHandler;
    class LinearMatrixHandler;
    class NonlinearMatrixHandler;
    class WcsMatrixHandler;
    class PsfMatrixHandler;

    boost::scoped_ptr<ModelImageHandler> _modelImageHandler;
    boost::scoped_ptr<LinearMatrixHandler> _linearMatrixHandler;
    boost::scoped_ptr<NonlinearMatrixHandler> _nonlinearMatrixHandler;
    boost::scoped_ptr<WcsMatrixHandler> _wcsMatrixHandler;
    boost::scoped_ptr<PsfMatrixHandler> _psfMatrixHandler;

};

class FourierModel::Factory  : public Model::Factory {
public:
    typedef boost::shared_ptr<Factory> Ptr;
    typedef boost::shared_ptr<const Factory> ConstPtr;
    
    virtual FourierModel * project(
        ParameterConstIterator linearParameterBegin,
        ParameterConstIterator const linearParameterEnd,
        ParameterConstIterator nonlinearParameterBegin,
        Kernel const & kernel,
        WcsConstPtr const & wcs,
        FootprintConstPtr const & footprint,
        double photFactor
    ) const;

    virtual int const getNonlinearParameterSize() const;

    virtual int const getMinLinearParameterSize() const {
        return _morphologyFactory->getMinLinearParameterSize();
    }

    virtual int const getMaxLinearParameterSize() const {
        return _morphologyFactory->getMaxLinearParameterSize();
    }

    Factory(
        Position::ConstPtr const & position, 
        Morphology::Factory::ConstPtr const & morphologyFactory
    ) : _position(position),
        _morphologyFactory(morphologyFactory)
    {}

private:
    Position::ConstPtr _position;
    Morphology::Factory::ConstPtr _morphologyFactory;
};

inline lsst::afw::image::PointD FourierModel::getPsfPosition() const { 
    return (*_transform)(_position->getCenter()); 
}

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_FOURIERMODEL_H
