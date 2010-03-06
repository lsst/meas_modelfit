#ifndef LSST_MEAS_MULTIFIT_SERSIC_MORPHOLOGY_H
#define LSST_MEAS_MULTIFIT_SERSIC_MORPHOLOGY_H

#include "lsst/meas/multifit/components/Morphology.h"
#include "lsst/meas/multifit/SersicCache.h"

namespace lsst {
namespace meas {
namespace multifit {
namespace components {

/**
 * Derived Morphology component for fitting small galaxies
 */
class SersicMorphology : public Morphology {
public:
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
    static SersicMorphology::Ptr createTemplate() {
        return SersicMorphology::Ptr(new SersicMorphology());
    }

    double const & getSersicIndex() const {
        return *(_getNonlinearParameterIter() + 3);
    }
    Cache::Functor::ConstPtr const & getSersicIndexFunctor() const {
        return _indexFunctor;
    }
protected:
    // Morphology *************************************************************
    virtual Morphology::Ptr create(
        boost::shared_ptr<ParameterVector const> const & linearParameters,
        ParameterConstIterator nonlinearParameterIter
    ) const;

    virtual int const getMinLinearParameterSize() const {return 1;}
    virtual int const getMaxLinearParameterSize() const {return 1;}
    virtual int const getNonlinearParameterSize() const {return 4;}

    virtual void _handleNonlinearParameterChange();

    /**
     * Default-constructor
     *
     * Used to create a template SersicMorphology
     */
    SersicMorphology() : Morphology() {}

    /**
     * Construct a SersicMorphology
     *
     * For use within ComponentModel
     * @sa Morphology::create()
     */
    SersicMorphology(
        boost::shared_ptr<ParameterVector const> const & linearParameters,
        ParameterConstIterator nonlinearParameterIter
    ) : Morphology(linearParameters, nonlinearParameterIter) {}

    Cache::Functor::ConstPtr _indexFunctor;
};


}}}} //namespace lsst::meas::multifit
#endif
