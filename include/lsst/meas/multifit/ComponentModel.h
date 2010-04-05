// -*- lsst-c++ -*-
/**
 * @file
 * Support for Models which are defined by an astrometry and a morphology
 * component
 */
#ifndef LSST_MEAS_MULTIFIT_COMPONENT_MODEL_H
#define LSST_MEAS_MULTIFIT_COMPONENT_MODEL_H

#include <boost/make_shared.hpp>

#include "lsst/afw/geom/Box.h"

#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/ModelProjection.h"
#include "lsst/meas/multifit/components/Astrometry.h"
#include "lsst/meas/multifit/components/FixedAstrometry.h"
#include "lsst/meas/multifit/components/Morphology.h"

namespace lsst {
namespace meas {
namespace multifit {

class ComponentModelProjection;

/**
 *  Derived Model defined by its a Astrometry and a Morphology components.
 *
 *  The number of nonlinear parameters is the sum of the Astrometry parameters
 *  and the Morphology parameters. When calling setNonlinearParameters on a
 *  ComponentModel, the input ParameterVector must be organized as 
 *  @sa Model
 *  @sa ComponentModelFactory
 *  @sa ComponentModelProjection
 */
class ComponentModel : public Model {
public:

    typedef boost::shared_ptr<ComponentModel> Ptr;
    typedef boost::shared_ptr<ComponentModel const> ConstPtr;

    virtual Footprint::Ptr computeProjectionFootprint(
        PsfConstPtr const & psf,
        WcsConstPtr const & wcs
    ) const;

    virtual lsst::afw::geom::BoxD computeProjectionEnvelope(
        PsfConstPtr const & psf,
        WcsConstPtr const & wcs
    ) const;

    virtual lsst::afw::geom::ellipses::Ellipse::Ptr computeBoundingEllipse() const;

    virtual Model::Ptr clone() const { 
        return Model::Ptr(new ComponentModel(*this)); 
    }

    /**
     *  Immutable access to the Astrometry component.
     */    
    components::Astrometry::ConstPtr getAstrometry() const {
        return _astrometry; 
    }

    /**
     *  Immutable access to the the Morphology component.
     */
    components::Morphology::ConstPtr getMorphology() const { 
        return _morphology; 
    }
    
    ComponentModel(
        Astrometry::ConstPtr const & astrometry, 
        Morphology::ConstPtr const & morphology,
    );
protected:

    virtual boost::shared_ptr<ModelProjection> makeProjection(
        PsfConstPtr const & psf,
        WcsConstPtr const & wcs,
        FootprintConstPtr const & footprint
    ) const;

private:
    friend class ComponentModelFactory;

    explicit ComponentModel(
        int linearParameterSize,
        components::Astrometry::ConstPtr const & astrometryTemplate,
        components::Morphology::ConstPtr const & morphologyTemplate
    );
    ComponentModel(ComponentModel const & model);
    void _initializeComponents(
        components::Astrometry::ConstPtr const & astrometryTemplate,
        components::Morphology::ConstPtr const & morphologyTemplate
    );


    /**
     * Iterator at the start of the nonlinear parameters relevant to the
     * astrometry component
     */
    ParameterConstIterator _getAstrometryParameterIter() const {
        return getNonlinearParameterIter();
    }
    /**
     * Iterator at the start of the nonlinear parameters relevant to the
     * morphology component
     */
    ParameterConstIterator _getMorphologyParameterIter() const {
        return getNonlinearParameterIter() + _astrometry->getParameterSize();
    }


    virtual void _handleLinearParameterChange();

    virtual void _handleNonlinearParameterChange();

    components::Astrometry::Ptr _astrometry;
    components::Morphology::Ptr _morphology;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_COMPONENT_MODEL_H
