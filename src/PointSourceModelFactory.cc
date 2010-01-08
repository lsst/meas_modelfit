#include "lsst/meas/multifit/PointSourceModelFactory.h"
#include <boost/make_shared.hpp>
namespace multifit = lsst::meas::multifit;

bool multifit::PointSourceModelFactory::declareMe() {
    static bool isDeclared = false;

    if(!isDeclared) {
        ModelFactory::declare(
            getName(), 
            boost::make_shared<PointSourceModelFactory const>()
        );
        isDeclared = true;
    }

    return isDeclared;
}

multifit::Model::Ptr multifit::PointSourceModelFactory::makeModel(
    double flux, lsst::afw::geom::Point2D const & position
) {
    ParameterVector linear(1);
    linear << flux;
    return ComponentModelFactory::makeModel(
        linear, 
        position.asVector()
    ); 
}

namespace {
 bool isDeclared = multifit::PointSourceModelFactory::declareMe();
}
