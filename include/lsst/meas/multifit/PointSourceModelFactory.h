#ifndef LSST_MEAS_MULTIFIT_POINT_SOURCE_MODEL_FACTORY_H
#define LSST_MEAS_MULTIFIT_POINT_SOURCE_MODEL_FACTORY_H
#include <string>

#include "lsst/meas/multifit/ComponentModelFactory.h"
#include "lsst/meas/multifit/components/PointSourceMorphology.h"

namespace lsst {
namespace meas {
namespace multifit {

class PointSourceModelFactory : public ComponentModelFactory { 
public:
    PointSourceModelFactory() :
        ComponentModelFactory(
            boost::make_shared<components::Astrometry>(), 
            components::PointSourceMorphology::createTemplate()
        )
    {}

    static bool declareMe();
    Model::Ptr makeModel(double, lsst::afw::geom::Point2D const &);
private:
    static std::string getName() {
        return "PointSource";
    }
};

}}} //end namespace lsst::meas::multifit

#endif //LSST_MEAS_MULTIFIT_POINT_SOURCE_MODEL_FACTORY_H
