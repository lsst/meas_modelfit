// -*- lsst-c++ -*-
/**
 * @file
 *
 * Implementation of a specialized ModelFactory used to create static point 
 * source models
 */
#ifndef LSST_MEAS_MULTIFIT_POINT_SOURCE_MODEL_FACTORY_H
#define LSST_MEAS_MULTIFIT_POINT_SOURCE_MODEL_FACTORY_H
#include <string>

#include "lsst/meas/multifit/ComponentModelFactory.h"
#include "lsst/meas/multifit/components/PointSourceMorphology.h"

namespace lsst {
namespace meas {
namespace multifit {


/**
 * A factory for a PointSourceModel
 *
 * There is no class implementation for a PointSourceModel, as it is just a
 * special case of a ComponentModel. However, a factory
 * implementation for point-sources allows implementing a specialized
 * constructor, which is convenient for testing and usage.
 *
 * PointSourceModelFactory is a very lightweight wrapper over
 * ComponentModelFactory, initializing its base class using the static 
 * Astrometry component, and the PointSourceMorphology component.
 *
 * @sa lsst::meas::multifit::ComponentModelFactory
 * @sa lsst::meas::multifit::components::PointSourceMorphology
 * @sa lsst::meas::multifit::components::Astrometry
 */
class PointSourceModelFactory : public ComponentModelFactory { 
public:

    /**
     * Default construct a PointSourceModelFactory
     *
     * Initializes its base class with the static Astrometry component, and 
     * the PointSourceMorphology component
     */
    PointSourceModelFactory() :
        ComponentModelFactory(
            boost::make_shared<components::Astrometry>(), 
            components::PointSourceMorphology::createTemplate()
        )
    {}


    static bool declareMe();
    
    Model::Ptr makeModel(double, lsst::afw::geom::Point2D const &);
private:

    /**
     * Name used to register this factory in the ModelFactory registry
     */
    static std::string getName() {
        return "PointSource";
    }
};

}}} //end namespace lsst::meas::multifit

#endif //LSST_MEAS_MULTIFIT_POINT_SOURCE_MODEL_FACTORY_H
