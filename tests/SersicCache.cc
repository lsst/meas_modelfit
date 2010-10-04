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
#define BOOST_TEST_MODULE SersicCache

#include "boost/make_shared.hpp"
#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"

#include "ndarray.hpp"

#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/SersicCache.h"

using namespace std;
namespace multifit = lsst::meas::multifit;

BOOST_AUTO_TEST_CASE(SersicCache) {
    multifit::Cache::ConstPtr cache;
    try {
        cache = multifit::Cache::load("testCache");
    } catch (...) {
        lsst::pex::policy::DefaultPolicyFile file("meas_multifit", "SersicCache.paf", "tests");
        lsst::pex::policy::Policy pol;
        file.load(pol);
        cache = multifit::makeSersicCache(pol);
    }

    BOOST_CHECK(cache);

    double slowMin = cache->getSlowMin(), slowMax = cache->getSlowMax();
    double fastMin = cache->getFastMin(), fastMax = cache->getFastMax(); 
    multifit::Cache::Interpolator::ConstPtr functor;
    BOOST_CHECK_NO_THROW(functor = cache->getInterpolator(slowMin));
    BOOST_CHECK(functor);
    BOOST_CHECK_NO_THROW((*functor)(fastMin));
    BOOST_CHECK_NO_THROW((*functor)(fastMax-1.0));

    BOOST_CHECK_NO_THROW(functor = cache->getDerivativeInterpolator(slowMin));
    BOOST_CHECK(functor);
    BOOST_CHECK_NO_THROW((*functor)(fastMin));
    BOOST_CHECK_NO_THROW((*functor)(fastMax-1.0));

    BOOST_CHECK_NO_THROW(functor = cache->getInterpolator(slowMax -1.0));
    BOOST_CHECK(functor);
    BOOST_CHECK_NO_THROW((*functor)(fastMin));
    BOOST_CHECK_NO_THROW((*functor)(fastMax - 1.0));
    
    BOOST_CHECK_NO_THROW(functor = cache->getDerivativeInterpolator(slowMax- 1.0));
    BOOST_CHECK(functor);
    BOOST_CHECK_NO_THROW((*functor)(fastMin));
    BOOST_CHECK_NO_THROW((*functor)(fastMax - 1.0));

    double const * data = cache->getDataPoints().data();
    for(int i = 0; i < cache->getDataPoints().size(); ++i, ++data) {
        //check that there are no NaNs in the cache
        BOOST_CHECK_EQUAL(*data, *data);
    }
}
