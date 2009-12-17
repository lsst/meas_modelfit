#include "lsst/meas/multifit/PointSourceModelFactory.h"
#include <boost/make_shared.hpp>
namespace multifit = lsst::meas::multifit;

bool multifit::PointSourceModelFactory::registerMe() {
    static bool isRegistered = false;

    if(!isRegistered) {
        ModelFactory::registerFactory(
            getName(), 
            boost::make_shared<PointSourceModelFactory const>()
        );
        isRegistered = true;
    }

    return isRegistered;
}

multifit::Model::Ptr multifit::PointSourceModelFactory::makeModel(
    double flux, lsst::afw::geom::Point2D const & position
) {
    ParameterVector linear(1);
    linear << flux;
    return ComponentModelFactory::makeModel(
        1, linear.data(), position.asVector().data()
    ); 
}

namespace {
 bool isRegistered = multifit::PointSourceModelFactory::registerMe();
}
