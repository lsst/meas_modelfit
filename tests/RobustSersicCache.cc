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
#define BOOST_TEST_MODULE RobustSersicCache

#include "boost/make_shared.hpp"
#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"

#include "ndarray.hpp"

#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/RobustSersicCache.h"

#include "lsst/pex/logging/Log.h"

using namespace std;
namespace multifit = lsst::meas::multifit;
namespace pexLog = lsst::pex::logging;

BOOST_AUTO_TEST_CASE(RobustSersicCache) {
    lsst::pex::logging::Log::getDefaultLog().setThresholdFor("lsst.meas.multifit", -10);

    multifit::Cache::ConstPtr cache;
    try {
        cache = multifit::Cache::load("testRobustCache");
    } catch (...) {
        lsst::pex::policy::Policy pol;
        cache = multifit::makeRobustSersicCache(pol);
    }

    BOOST_CHECK(cache);

    lsst::afw::geom::BoxD bounds = cache->getParameterBounds();
    
    multifit::Cache::Functor::ConstPtr functor;
    BOOST_CHECK_NO_THROW(functor = cache->getRowFunctor(bounds.getMinY()));
    BOOST_CHECK(functor);
    BOOST_CHECK_NO_THROW((*functor)(bounds.getMinX()));
    BOOST_CHECK_NO_THROW((*functor)(bounds.getMaxX()-0.1));

    BOOST_CHECK_NO_THROW(functor = cache->getRowFunctor(bounds.getMaxY() - 0.1));
    BOOST_CHECK(functor);
    BOOST_CHECK_NO_THROW((*functor)(bounds.getMinX()));
    BOOST_CHECK_NO_THROW((*functor)(bounds.getMaxX()-0.1));

    BOOST_CHECK_NO_THROW(functor = cache->getColFunctor(bounds.getMinX()));
    BOOST_CHECK(functor);
    BOOST_CHECK_NO_THROW((*functor)(bounds.getMinY()));
    BOOST_CHECK_NO_THROW((*functor)(bounds.getMaxY()-0.1));
    
    BOOST_CHECK_NO_THROW(functor = cache->getColFunctor(bounds.getMaxX() - 0.1));
    BOOST_CHECK(functor);
    BOOST_CHECK_NO_THROW((*functor)(bounds.getMinY()));
    BOOST_CHECK_NO_THROW((*functor)(bounds.getMaxY()-0.1));

    double const * data = cache->getDataPoints().data();
    for(int i = 0; i < cache->getDataPoints().size(); ++i, ++data) {
        //check that there are no NaNs in the cache
        BOOST_CHECK_EQUAL(*data, *data);
    }
}
