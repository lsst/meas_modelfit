#ifndef LSST_MEAS_MULTIFIT_COMPONENT_MODEL_H
#define LSST_MEAS_MULTIFIT_COMPONENT_MODEL_H

#include <boost/make_shared.hpp>

#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/ModelProjection.h"
#include "lsst/meas/multifit/components/Astrometry.h"
#include "lsst/meas/multifit/components/Morphology.h"

namespace lsst {
namespace meas {
namespace multifit {

class ComponentModelProjection;

/**
 *  \brief A Model defined by a Astrometry object and a Morphology.
 *
 *  \sa Model
 *  \sa ComponentModelFactory
 *  \sa ComponentModelProjection
 */
class ComponentModel : public Model {
public:

    typedef boost::shared_ptr<ComponentModel> Ptr;
    typedef boost::shared_ptr<ComponentModel const> ConstPtr;

    /**
     *  \brief Create a Footprint that would contain a projection of the Model.
     */
    virtual Footprint::Ptr computeProjectionFootprint(
        Kernel const & kernel,
        Wcs::ConstPtr const & wcs,
        double photFactor
    ) const;

    /**
     *  \brief Create an image-coordinate bounding box that would contain a projection of the Model.
     */
    virtual agl::BoxD computeProjectionEnvelope(
        Kernel const & kernel,
        Wcs::ConstPtr const & wcs,
        double photFactor
    ) const;

    /**
     *  \brief Create an ra/dec bounding ellipse for the Model (does not include PSF broadening).
     */
    virtual lsst::afw::math::ellipses::Ellipse::Ptr computeBoundingEllipse() const;

    /**
     *  \brief Create a new ComponentModel with the same type and parameters.
     *
     *  Associated ModelProjection objects are not copied or shared;
     *  the new Model will not have any associated ModelProjections.
     */
    virtual Model::Ptr clone() const { return Model::Ptr(new ComponentModel(*this)); }

    /**
     *  \brief Return the Astrometry component of the Model.
     */
    components::Astrometry::ConstPtr getAstrometry() const { return _astrometry; }

    /**
     *  \brief Return the Morphology component of the Model.
     */
    components::Morphology::ConstPtr getMorphology() const { return _morphology; }

protected:

    /** 
     *  \brief Create a ModelProjection object associated with this.
     */
    virtual boost::shared_ptr<ModelProjection> makeProjection(
        Kernel const & kernel,
        Wcs::ConstPtr const & wcs,
        Footprint::ConstPtr const & footprint,
        double photFactor,
        int activeProducts
    ) const;

private:

    friend class ComponentModelFactory;

    /// \brief Initialize the ComponentModel (does not set parameters).
    explicit ComponentModel(
        int linearParameterSize,
        components::Astrometry::ConstPtr const & astrometryTemplate,
        components::Morphology::ConstPtr const & morphologyTemplate
    );

    /// \brief Deep-copy the Model.
    ComponentModel(ComponentModel const & model);

    /// \brief Shared initialization code for ComponentModel constructors.
    void _construct(
        components::Astrometry::ConstPtr const & astrometryTemplate,
        components::Morphology::ConstPtr const & morphologyTemplate
    );

    /// \brief Notifies the Morphology object that the linear parameters have changed.
    virtual void _handleLinearParameterChange();

    /// \brief Notifies the Astrometry and Morphology objects that the nonlinear parameters have changed.
    virtual void _handleNonlinearParameterChange();

    components::Astrometry::Ptr _astrometry;
    components::Morphology::Ptr _morphology;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_COMPONENT_MODEL_H
