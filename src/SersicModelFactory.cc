// -*- lsst-c++ -*-
/**
 * @file
 * Implementation of SersicModelFactory
 */

#include "lsst/meas/multifit/SersicModelFactory.h"
#include <boost/make_shared.hpp>
namespace multifit = lsst::meas::multifit;
/**
 * Add this ModelFactory type to the ModelFactory registry
 *
 * The factory is registered under the name 'Sersic'
 *
 * @sa lsst::meas::multifit::ModelFactory::declare
 */
bool multifit::SersicModelFactory::declareMe() {
    static bool isDeclared = false;

    if(!isDeclared) {
        ModelFactory::declare(
            getName(), 
            boost::make_shared<SersicModelFactory const>()
        );
        isDeclared = true;
    }

    return isDeclared;
}

/**
 * Specialized Model constructor
 *
 * @param flux
 * @param position model centroid
 * @param ellipseCore an parametrization of an ellipse, will be converted to
 *     LogShear parametrization.
 * @param sersicIndex
 */
multifit::Model::Ptr multifit::SersicModelFactory::makeModel(
    double const & flux, 
    lsst::afw::geom::Point2D const & position,
    lsst::afw::geom::ellipses::Core::ConstPtr const & ellipseCore,
    double const & sersicIndex
) {
    ParameterVector linear(1), nonlinear(6);
    linear << flux;
    nonlinear << position.asVector(), 
        static_cast<lsst::afw::geom::ellipses::LogShear>(*ellipseCore).getVector(),
        sersicIndex;

    return ComponentModelFactory::makeModel(
        linear, 
        nonlinear
    ); 
}

/**
 * Specialized Model constructor
 *
 * @param flux
 * @param ellipse an parametrization of an ellipse, will be converted to
 *     LogShear parametrization.
 * @param sersicIndex
 */
multifit::Model::Ptr multifit::SersicModelFactory::makeModel(
    double const & flux, 
    lsst::afw::geom::ellipses::Ellipse::ConstPtr const & ellipse,
    double const & sersicIndex
) {
    ParameterVector linear(1), nonlinear(6);
    linear << flux;
    nonlinear << 
        static_cast<lsst::afw::geom::ellipses::LogShearEllipse>(*ellipse).getVector(),
        sersicIndex;

    return ComponentModelFactory::makeModel(
        linear, 
        nonlinear
    ); 
}

namespace {
    
    bool isDeclared = multifit::SersicModelFactory::declareMe();
}
