#ifndef LSST_MEAS_MULTIFIT_COMPONENTS_SERSIC_MORPHOLOGY_PROJECTION_H
#define LSST_MEAS_MULTIFIT_COMPONENTS_SERSIC_MORPHOLOGY_PROJECTION_H

#include "lsst/meas/multifit/components/FourierMorphologyProjection.h"
#include "lsst/meas/multifit/components/SersicMorphology.h"
#include "lsst/meas/multifit/EllipseGridTransform.h"
#include <ndarray/fft.hpp>
namespace lsst {
namespace meas {
namespace multifit {
namespace components {

class SersicMorphologyProjection : public FourierMorphologyProjection {
public:
    typedef boost::shared_ptr<SersicMorphologyProjection> Ptr;
    typedef boost::shared_ptr<SersicMorphologyProjection const> ConstPtr;
    
    virtual lsst::afw::geom::Extent2I getDimensions() const;

    virtual ParameterJacobianMatrixPtr computeProjectedParameterJacobian() const;
    virtual TransformJacobianMatrixPtr computeTransformParameterJacobian() const; 

    virtual ndarray::FourierArray<Pixel,3,3> computeLinearParameterDerivative();
    virtual ndarray::FourierArray<Pixel,3,3> computeProjectedParameterDerivative();

    SersicMorphology::ConstPtr getMorphology() const {
        return boost::static_pointer_cast<SersicMorphology const>(
            MorphologyProjection::getMorphology()
        );
    }
protected:
    friend class SersicMorphology;
    typedef ndarray::FourierTransform<Pixel, 2> FFT;

    virtual void _handleLinearParameterChange();
    virtual void _handleNonlinearParameterChange();

    SersicMorphologyProjection(
        SersicMorphology::ConstPtr const & morphology,
        lsst::afw::geom::Extent2I const & kernelDimensions,
        lsst::afw::geom::AffineTransform::ConstPtr const & transform
    ) : FourierMorphologyProjection(morphology, kernelDimensions, transform), 
        _dimensions(),
        _validProducts(0)
    {
        _recomputeDimensions();
    }

    EllipseGridTransform::ConstPtr computeEllipseGridTransform() const;
    ndarray::FourierArray<Pixel, 3, 3> _projectedParameterDerivative;
    ndarray::FourierArray<Pixel, 3, 3> _linearParameterDerivative;

private:
    void _recomputeDimensions();
    lsst::afw::geom::Extent2I _dimensions;

    int _validProducts;
    enum Products {
        LINEAR_PARAMETER_DERIVATIVE=1,
        PROJECTED_PARAMETER_DERIVATIVE=2,
    };
};

}}}}
#endif
