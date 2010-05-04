// -*- LSST-C++ -*-
#include <iostream>
#include <cmath>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE FixedAstrometry

#include "boost/make_shared.hpp"
#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"

#include "ndarray.hpp"
#include "lsst/afw/geom.h"
#include "lsst/meas/multifit/ComponentModel.h"
#include "lsst/meas/multifit/components/Astrometry.h"
#include "lsst/meas/multifit/components/FixedAstrometry.h"
#include "lsst/meas/multifit/components/PointSourceMorphology.h"

using namespace std;
namespace multifit = lsst::meas::multifit;
namespace components = multifit::components;
namespace geom = lsst::afw::geom;

BOOST_AUTO_TEST_CASE(FixedAstrometry) {
    geom::Point2D point = lsst::afw::geom::makePointD(12, 12);
    components::Astrometry::Ptr astrometry(
        new components::Astrometry(point)
    );
    components::Astrometry::Ptr fixedAstrometry(
        new components::FixedAstrometry(point)
    );

    BOOST_CHECK_EQUAL(astrometry->computePosition(), point);
    BOOST_CHECK_EQUAL(fixedAstrometry->computePosition(), point);
    
    double flux = 1.0;
    components::Morphology::Ptr morphology = 
        components::PointSourceMorphology::create(flux);

    multifit::ComponentModel::Ptr model;

    model = multifit::ComponentModel::create(astrometry, morphology);
    BOOST_CHECK_EQUAL(model->getNonlinearParameterSize(), 2);

    model = multifit::ComponentModel::create(fixedAstrometry, morphology);
    BOOST_CHECK_EQUAL(model->getNonlinearParameterSize(), 0);
}
