// -*- lsst-c++ -*-
/**
 * @file
 * Declaration of class PointSourceMorphology
 */
#ifndef LSST_MEAS_MULTIFIT_COMPONENTS_POINT_SOURCE_MORPHOLOGY_H
#define LSST_MEAS_MULTIFIT_COMPONENTS_POINT_SOURCE_MORPHOLOGY_H

#include <boost/make_shared.hpp>

#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/components/Morphology.h"

namespace lsst {
namespace meas {
namespace multifit {
namespace components {

/**
 * Derived Morphology component for fitting static point-sources
 */
class PointSourceMorphology : public Morphology {
public:
    typedef boost::shared_ptr<PointSourceMorphology> Ptr;
    typedef boost::shared_ptr<PointSourceMorphology const> ConstPtr;


    virtual lsst::afw::geom::ellipses::Core::Ptr computeBoundingEllipseCore() const {
        return boost::make_shared<lsst::afw::geom::ellipses::LogShear>();
    }

    virtual MorphologyProjection::Ptr makeProjection(
        lsst::afw::geom::Extent2I const & kernelDimensions,
        lsst::afw::geom::AffineTransform::ConstPtr const & transform
    ) const;

    /**
     * Named PointSourceMorphology constructor    
     */
    static PointSourceMorphology::Ptr createTemplate() {
        return PointSourceMorphology::Ptr(new PointSourceMorphology());
    }       

protected:

    /**
     *  Default-construct a Morphology object to be used as a template.
     */
    PointSourceMorphology() : Morphology() {}

    /**
     *  Construct a Morphology object for use inside a ComponentModel.
     *
     *  @sa Morphology::create()
     */
    PointSourceMorphology(
        boost::shared_ptr<ParameterVector const> const & linearParameterVector,
        ParameterConstIterator nonlinearParameterIter
    ) : Morphology(linearParameterVector,nonlinearParameterIter) {}

    virtual Morphology::Ptr create(
        boost::shared_ptr<ParameterVector const> const & linearParameterVector,
        ParameterConstIterator nonlinearParameterIter 
    ) const;

    //derived Morphology template mode functions
    virtual int const getMinLinearParameterSize() const { return 1; }
    virtual int const getMaxLinearParameterSize() const { return 1; }
    virtual int const getNonlinearParameterSize() const { return 0; }

};

}}}} // namespace lsst::meas::multifit::components

#endif // !LSST_MEAS_MULTIFIT_COMPONENTS_POINT_SOURCE_MORPHOLOGY_H
