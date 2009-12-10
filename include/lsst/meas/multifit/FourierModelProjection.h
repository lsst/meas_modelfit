#ifndef LSST_MEAS_MULTIFIT_FOURIER_MODEL_PROEJECTION_H
#define LSST_MEAS_MULTIFIT_FOURIER_MODEL_PROEJECTION_H

#include <ndarray/fft_fwd.hpp>

#include "lsst/afw/math/AffineTransform.h"
#include "lsst/afw/math/ConvolutionVisitor.h"
#include "lsst/afw/image/Utils.h"

#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/ComponentModelProjection.h"
#include "lsst/meas/multifit/components/FourierMorphologyProjection.h"

namespace lsst {
namespace meas {
namespace multifit {

class FourierModelProjection : public ComponentModelProjection {
public:

    typedef boost::shared_ptr<FourierModelProjection> Ptr;
    typedef boost::shared_ptr<FourierModelProjection const> ConstPtr;

    /// \brief Return the number of parameters that specify the PSF.
    virtual int const getPsfParameterSize() const;

    /// \brief Return the Model instance this is a projection of.
    ComponentModel::ConstPtr getModel() const {
        return boost::static_pointer_cast<ComponentModel const>(
            ModelProjection::getModel()
        );
    }

    /// \brief Return the MorphologyProjection object this projection is based on.
    components::FourierMorphologyProjection::ConstPtr getMorphologyProjection() const { 
        return boost::static_pointer_cast<components::FourierMorphologyProjection const>(
            getMorphologyProjection()
        );
    }

    ~FourierModelProjection();

protected:

    friend class multifit::ComponentModel;

    /**
     *  @name ConvolvableImplementation
     */
    //@{
    virtual void convolve(lsst::afw::math::FourierConvolutionVisitor const & visitor);

    virtual bool isConvolved() const { return _kernelVisitor; }
    //@}

    /**
     *  @name ProtectedProductComputers
     */
    //@{
    virtual void _computeLinearParameterDerivative(ndarray::Array<Pixel,2,2> const & matrix);
    virtual void _computePsfParameterDerivative(ndarray::Array<Pixel,2,2> const & matrix);
    virtual void _computeTranslationDerivative(ndarray::Array<Pixel,2,2> const & matrix);
    virtual void _computeProjectedParameterDerivative(ndarray::Array<Pixel,2,2> const & matrix);
    //@}

    virtual int _enableProducts(int toAdd);
    virtual int _disableProducts(int toRemove);

    virtual void _handleLinearParameterChange();
    virtual void _handleNonlinearParameterChange();

    components::FourierMorphologyProjection::Ptr _getMorphologyProjection() { 
        return boost::static_pointer_cast<components::FourierMorphologyProjection>(
            _getMorphologyProjection()
        );
    }

private:

    typedef ndarray::FourierTransform<Pixel,2> FFT;
    
    FourierModelProjection(
        ComponentModel::ConstPtr const & model,
        Kernel const & kernel,
        Wcs::ConstPtr const & wcs,
        Footprint::ConstPtr const & footprint,
        double photFactor,
        int activeProducts = 0
    );

    void _setDimensions();

    void _applyKernel(
        ndarray::FourierArray<Pixel,3,3>::Iterator iter, 
        ndarray::FourierArray<Pixel,3,3>::Iterator const & end,
        ndarray::FourierArray<Pixel const,2,2> const & kernel
    ) const;

    void _applyKernel(
        ndarray::FourierArray<Pixel,3,3>::Iterator iter, 
        ndarray::FourierArray<Pixel,3,3>::Iterator const & end
    ) const {
        _applyKernel(iter,end,_kernelVisitor->getFourierImage());
    }

    lsst::afw::math::FourierConvolutionVisitor::Ptr _kernelVisitor; ///< PSF/Kernel information.
    WindowedFootprint::Ptr _wf; ///< maps footprint to handler output arrays
    lsst::afw::image::BBox _outerBBox; ///< bounding box of padded arrays relative to exposure
    lsst::afw::image::BBox _innerBBox; ///< bounding box of unpadded arrays relative to _outerBBox

    class Shifter;
    class LinearMatrixHandler;
    class NonlinearMatrixHandler;
    class PsfMatrixHandler;

    boost::scoped_ptr<Shifter> _shifter;
    boost::scoped_ptr<LinearMatrixHandler> _linearMatrixHandler;
    boost::scoped_ptr<NonlinearMatrixHandler> _nonlinearMatrixHandler;
    boost::scoped_ptr<PsfMatrixHandler> _psfMatrixHandler;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_FOURIER_MODEL_PROEJECTION_H

