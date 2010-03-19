// -*- LSST-C++ -*-
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
    multifit::SersicCache::ConstPtr cache = multifit::SersicCache::getInstance();

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
