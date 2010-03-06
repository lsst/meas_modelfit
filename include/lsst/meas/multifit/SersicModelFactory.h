// -*- lsst-c++ -*-
/**
 * @file
 *
 * Implementation of a specialized ModelFactory used to create sersic models
 */
#ifndef LSST_MEAS_MULTIFIT_SERSIC_MODEL_FACTORY_H
#define LSST_MEAS_MULTIFIT_SERSIC_MODEL_FACTORY_H
#include <string>

#include "lsst/meas/multifit/ComponentModelFactory.h"
#include "lsst/meas/multifit/components/SersicMorphology.h"

namespace lsst {
namespace meas {
namespace multifit {


/**
 * A factory for a SersicModel
 *
 * There is no class implementation for a SersicModel, as it is just a
 * special case of a ComponentModel. However, a factory
 * specialization allows implementing a specialized
 * constructor, which is convenient for testing and usage.
 *
 * SersicModelFactory is a very lightweight wrapper over
 * ComponentModelFactory, initializing its base class using the static 
 * Astrometry component, and the PointSourceMorphology component.
 *
 * @sa lsst::meas::multifit::ComponentModelFactory
 * @sa lsst::meas::multifit::components::SeriscMorphology
 * @sa lsst::meas::multifit::components::Astrometry
 */
class SersicModelFactory : public ComponentModelFactory { 
public:

    /**
     * Default construct a PointSourceModelFactory
     *
     * Initializes its base class with the static Astrometry component, and 
     * the PointSourceMorphology component
     */
    SersicModelFactory() :
        ComponentModelFactory(
            boost::make_shared<components::Astrometry>(), 
            components::SersicMorphology::createTemplate()
        )
    {}


    static bool declareMe();
    
    Model::Ptr makeModel(
        double const & flux, 
        lsst::afw::geom::Point2D const & centroid, 
        lsst::afw::geom::ellipses::Core::ConstPtr const & ellipseCore,
        double const & sesricIndex
    );
    Model::Ptr makeModel(
        double const & flux,
        lsst::afw::geom::ellipses::Ellipse::ConstPtr const & ellipse,
        double const & sersicIndex
    );
private:

    /**
     * Name used to register this factory in the ModelFactory registry
     */
    static std::string getName() {
        return "Sersic";
    }
};

}}} //end namespace lsst::meas::multifit

#endif //LSST_MEAS_MULTIFIT_POINT_SOURCE_MODEL_FACTORY_H
