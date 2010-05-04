#ifndef LSST_MEAS_MULTIFIT_SERSIC_MORPHOLOGY_H
#define LSST_MEAS_MULTIFIT_SERSIC_MORPHOLOGY_H

#include "lsst/meas/multifit/components/Morphology.h"
#include "lsst/meas/multifit/SersicCache.h"

namespace lsst {
namespace meas {
namespace multifit {
namespace components {

class SersicMorphologyProjection;

/**
 * Derived Morphology component for fitting small galaxies
 */
class SersicMorphology : public Morphology {
public:
    enum LinearParameters {
        FLUX=0,
        LINEAR_SIZE
    };
    enum NonlinearParameters {
        GAMMA1=0,
        GAMMA2,
        KAPPA,
        SERSIC_INDEX,
        NONLINEAR_SIZE
    };

    typedef boost::shared_ptr<SersicMorphology> Ptr;
    typedef boost::shared_ptr<SersicMorphology const> ConstPtr;

    virtual lsst::afw::geom::ellipses::BaseCore::Ptr computeBoundingEllipseCore() const;

    virtual MorphologyProjection::Ptr makeProjection(
        lsst::afw::geom::Extent2I const & kernelSize,
        lsst::afw::geom::AffineTransform::ConstPtr const & transform
    ) const;
  
    /**
     * Named SersicMorphology constructor
     */
    static SersicMorphology::Ptr create(
        Parameter const & flux,
        lsst::afw::geom::ellipses::BaseCore const & ellipse, 
        Parameter const & sersicIndex
    ) { 
        lsst::afw::geom::ellipses::LogShear logShear(ellipse);
        ParameterVector linear(LINEAR_SIZE), nonlinear(NONLINEAR_SIZE);
        linear << flux;
        nonlinear << logShear.getVector(), sersicIndex;
        return SersicMorphology::Ptr(            
            new SersicMorphology(
                boost::make_shared<ParameterVector const>(linear),
                boost::make_shared<ParameterVector const>(nonlinear)
            )
        );
    }

    Parameter const & getFlux() const {
        return *(_getLinearParameters()->data() + FLUX);
    }
    Parameter const & getSersicIndex() const {
        return *(beginNonlinear() + SERSIC_INDEX);
    }

protected:
    friend class SersicMorphologyProjection;

    virtual Morphology::Ptr create(
        boost::shared_ptr<ParameterVector const> const & linearParameters,
        boost::shared_ptr<ParameterVector const> const & nonlinearParameters,
        size_t const & start=0
    ) const;

    virtual void _handleNonlinearParameterChange();
    virtual int const getNonlinearParameterSize() const {return NONLINEAR_SIZE;}


    multifit::Cache::Functor::ConstPtr const & getSersicIndexFunctor() const {
        return _indexFunctor;
    }

    /**
     * Construct a SersicMorphology
     *
     * For use within ComponentModel
     * @sa Morphology::create()
     */
    SersicMorphology(
        boost::shared_ptr<ParameterVector const> const & linearParameters,
        boost::shared_ptr<ParameterVector const> const & nonlinearParameters,
        size_t const & start =0
    ) : Morphology(linearParameters, nonlinearParameters, start) {
        _handleNonlinearParameterChange();
    }

    Cache::Functor::ConstPtr _indexFunctor;
};


}}}} //namespace lsst::meas::multifit
#endif
