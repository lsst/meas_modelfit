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
    multifit::SersicCache::ConstPtr cache;
    try {
        cache = multifit::SersicCache::load("testCache");
    } catch (...) {
        lsst::pex::policy::DefaultPolicyFile file("meas_multifit", "SersicCache.paf", "tests");
        lsst::pex::policy::Policy pol;
        file.load(pol);
        cache = multifit::SersicCache::make(pol);
    }

    BOOST_CHECK(cache);

    double const eps1 = 1E-8;
    double const eps2 = 1E-6;

    double sersicMin = cache->getSersicMin(), sersicMax = cache->getSersicMax();
    double kMin = cache->getKMin(), kMax = cache->getKMax(); 

    for (double n = sersicMin; n <= sersicMax; n += 0.125) {
        double q = cache->convertSersicToParameter(n);
        double dq_dn = cache->differentiateSersicToParameter(n);
        double dn_dq = cache->differentiateParameterToSersic(q);
        double dq_dn_numeric = (cache->convertSersicToParameter(n + eps1) 
                                - cache->convertSersicToParameter(n - eps1)) / (2.0 * eps1);
        double dn_dq_numeric = (cache->convertParameterToSersic(q + eps1) 
                                - cache->convertParameterToSersic(q - eps1)) / (2.0 * eps1);
        BOOST_CHECK(std::fabs(cache->convertParameterToSersic(q) - n) < eps1);
        if (n > sersicMin && n < sersicMax) {
            BOOST_CHECK(std::fabs(dq_dn - dq_dn_numeric) < eps2);
            BOOST_CHECK(std::fabs(dn_dq - dn_dq_numeric) < eps2);
        }
    }

    BOOST_CHECK(sersicMin <= cache->convertParameterToSersic(1E6));
    BOOST_CHECK(sersicMin <= cache->convertParameterToSersic(1E-6));

    BOOST_CHECK(sersicMax >= cache->convertParameterToSersic(1E6));
    BOOST_CHECK(sersicMax >= cache->convertParameterToSersic(1E-6));

    multifit::SersicCache::Interpolator::ConstPtr functor;
    BOOST_CHECK_NO_THROW(functor = cache->getInterpolator(sersicMin));
    BOOST_CHECK(functor);
    BOOST_CHECK_NO_THROW((*functor)(kMin));
    BOOST_CHECK_NO_THROW((*functor)(kMax-1.0));

    BOOST_CHECK_NO_THROW(functor = cache->getDerivativeInterpolator(sersicMin));
    BOOST_CHECK(functor);
    BOOST_CHECK_NO_THROW((*functor)(kMin));
    BOOST_CHECK_NO_THROW((*functor)(kMax-1.0));

    BOOST_CHECK_NO_THROW(functor = cache->getInterpolator(sersicMax -1.0));
    BOOST_CHECK(functor);
    BOOST_CHECK_NO_THROW((*functor)(kMin));
    BOOST_CHECK_NO_THROW((*functor)(kMax - 1.0));
    
    BOOST_CHECK_NO_THROW(functor = cache->getDerivativeInterpolator(sersicMax- 1.0));
    BOOST_CHECK(functor);
    BOOST_CHECK_NO_THROW((*functor)(kMin));
    BOOST_CHECK_NO_THROW((*functor)(kMax - 1.0));

    double const * data = cache->getDataPoints().data();
    for(int i = 0; i < cache->getDataPoints().size(); ++i, ++data) {
        //check that there are no NaNs in the cache
        BOOST_CHECK_EQUAL(*data, *data);
    }
}
