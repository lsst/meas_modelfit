#ifndef LSST_MEAS_MULTIFIT_COMPONENTS_FIXED_NONLINEAR_MORPHOLOGY
#define LSST_MEAS_MULTIFIT_COMPONENTS_FIXED_NONLINEAR_MORPHOLOGY

#include "lsst/meas/multifit/components/Morphology.h"

namespace lsst{
namespace meas {
namespace multifit {
namespace components {

class FixedNonlinearMorphology : public Morphology {
public:
    typedef boost::shared_ptr<FixedNonlinearMorphology> Ptr;
    typedef boost::shared_ptr<FixedNonlinearMorphology const> ConstPtr;

    static Ptr create(Morphology const & base) {
        boost::shared_ptr<ParameterVector> nonlinearParameters;
        if (base.getNonlinearParameterSize() != 0) {
            nonlinearParameters.reset(
                new ParameterVector(base.getNonlinearParameterSize())
            );
            std::copy(
                base.beginNonlinear(), 
                base.endNonlinear(), 
                nonlinearParameters->data()
            );            
        } else {
            nonlinearParameters.reset(new ParameterVector());
        }
        
        Ptr fixed(new FixedNonlinearMorphology(base, nonlinearParameters));
        return fixed;
    }
    virtual lsst::afw::geom::ellipses::Core::Ptr computeBoundingEllipseCore() const {
        return _base->computeBoundingEllipseCore();
    };

    virtual int const getNonlinearParameterSize() const {return 0;}
    
    virtual MorphologyProjection::Ptr makeProjection (
        lsst::afw::geom::Extent2I const & kernelSize,
        lsst::afw::geom::AffineTransform::ConstPtr const & transform
    ) const {
        return _base->makeProjection(kernelSize, transform);
    }

    virtual Morphology::Ptr create(
        boost::shared_ptr<ParameterVector const> const &,
        boost::shared_ptr<ParameterVector const> const &,
        size_t const &
    ) const {
        Ptr fixed(new FixedNonlinearMorphology(*_base, getNonlinearParameters()));
        return fixed;
    }
protected:    
    Morphology::Ptr _base;

    virtual void _handleLinearParameterChange() {_base->_handleLinearParameterChange();}
    virtual void _handleNonlinearParameterChange() {};

    FixedNonlinearMorphology(
        Morphology const & base,
        boost::shared_ptr<ParameterVector const> const & nonlinearParameters
    ) : Morphology(base.getLinearParameters(), nonlinearParameters),
        _base(base.create(base.getLinearParameters(), nonlinearParameters)) {}
};

}}}}
#endif
