#ifndef LSST_MEAS_MULTIFIT_COMPONENTS_POINT_SOURCE_MORPHOLOGY_H
#define LSST_MEAS_MULTIFIT_COMPONENTS_POINT_SOURCE_MORPHOLOGY_H

#include <boost/make_shared.hpp>

#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/components/Morphology.h"

namespace lsst {
namespace meas {
namespace multifit {
namespace components {

class PointSourceMorphology : public Morphology {
public:
    typedef boost::shared_ptr<PointSourceMorphology> Ptr;
    typedef boost::shared_ptr<PointSourceMorphology const> ConstPtr;

    // --- Template-mode functionality ----------------------------------------------------------------------

    /// \brief Return the minimum number of linear parameters.
    virtual int const getMinLinearParameterSize() const { return 1; }

    /// \brief Return the maximum number of linear parameters.
    virtual int const getMaxLinearParameterSize() const { return 1; }

    /// \brief Return the number of (nonlinear) morphology parameters.
    virtual int const getMorphologyParameterSize() const { return 0; }

    // --- In-model functionality ---------------------------------------------------------------------------

    /// \brief Return an ellipse core that bounds the morphology.
    virtual lsst::afw::geom::ellipses::Core::Ptr computeBoundingEllipseCore() const {
        return boost::make_shared<lsst::afw::geom::ellipses::LogShear>();
    }

    /**
     *  \brief Create a new MorphologyProjection object.
     *
     *  Typically used only by the owning ComponentModel.
     */
    virtual MorphologyProjection::Ptr makeProjection(
        lsst::afw::geom::Extent2I const & kernelDimensions,
        lsst::afw::geom::AffineTransform::ConstPtr const & transform
    ) const;

    static PointSourceMorphology::Ptr createTemplate() {
        return PointSourceMorphology::Ptr(new PointSourceMorphology());
    }       

protected:

    /**
     *  \brief Default-construct a Morphology object to be used as a template.
     */
    PointSourceMorphology() : Morphology() {}

    /**
     *  \brief Construct a Morphology object for use inside a ComponentModel.
     *
     *  \sa Morphology::create()
     */
    PointSourceMorphology(
        boost::shared_ptr<ParameterVector const> const & linearParameterVector,
        ParameterConstIterator morphologyParameterIter
    ) : Morphology(linearParameterVector,morphologyParameterIter) {}

    /**
     *  \brief Construct a new Morphology using this as a template.
     *
     *  Typically used only by ComponentModel.
     *
     *  \param linearParameterVector The owning ComponentModel's linear parameter vector
     *  \param morphologyParameterIter Iterator to the first Morphology-specific
     *   parameter in the owning ComponentModel's nonlinear parameter vector.
     */
    virtual Morphology::Ptr create(
        boost::shared_ptr<ParameterVector const> const & linearParameterVector,
        ParameterConstIterator morphologyParameterIter 
    ) const;

};

}}}} // namespace lsst::meas::multifit::components

#endif // !LSST_MEAS_MULTIFIT_COMPONENTS_POINT_SOURCE_MORPHOLOGY_H
