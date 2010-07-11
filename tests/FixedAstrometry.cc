// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
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
