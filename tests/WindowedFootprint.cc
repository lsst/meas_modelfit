// -*- LSST-C++ -*-
#include <iostream>
#include <cmath>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE WindowedFootprint

#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"

#include "ndarray.hpp"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/meas/multifit/WindowedFootprint.h"

namespace multifit = lsst::meas::multifit;
namespace detection = lsst::afw::detection;
namespace geom = lsst::afw::geom;
namespace image = lsst::afw::image;

BOOST_AUTO_TEST_CASE(WindowedFootprintBasic) {
    geom::PointI min = geom::PointI::makeXY(3,4);
    int size = 20;
    geom::PointI max(min + geom::ExtentI(size-1));

    detection::Footprint fp(
        image::BBox(image::PointI(min.getX(), min.getY()), size, size)
    );
    image::BBox imageBBox = fp.getBBox();

    //test trivial case:
    //make a windowed fp that includes the entire fp
    geom::BoxI geomBBox = geom::BoxI(min, max);

    multifit::WindowedFootprint windowedFp(fp, geomBBox);
    BOOST_CHECK(fp.getNpix() != 0);
    BOOST_CHECK_EQUAL(fp.getNpix(), windowedFp.getNpix());
    BOOST_CHECK_EQUAL(windowedFp.getWindow().getWidth(), fp.getBBox().getWidth());
    BOOST_CHECK_EQUAL(windowedFp.getWindow().getHeight(), fp.getBBox().getHeight());
    //shrink the window in a predicatble way
    //this is halving the size of the box, also the number of pixels
    geomBBox.grow(geom::ExtentI::makeXY(0, -size/4));

    windowedFp = multifit::WindowedFootprint(fp, geomBBox);
    BOOST_CHECK_EQUAL(windowedFp.getWindow().getWidth(), fp.getBBox().getWidth());
    BOOST_CHECK_EQUAL(windowedFp.getWindow().getHeight(), fp.getBBox().getHeight()/2);
    BOOST_CHECK_EQUAL(windowedFp.getNpix(), fp.getNpix()/2);
}

BOOST_AUTO_TEST_CASE(WindowFootprintCompress) {
    geom::PointI min = geom::PointI::makeXY(3,4);
    int size = 4;

    detection::Footprint fp(
        image::BBox(image::PointI(min.getX(), min.getY()), size, size)
    );
    image::BBox imageBBox = fp.getBBox();

    geom::BoxI geomBBox = geom::BoxI(min, geom::ExtentI::makeXY(size, size/2));
    multifit::WindowedFootprint windowedFp(fp, geomBBox);

    BOOST_WARN(windowedFp.getNpix()>0);
    ndarray::Array<double, 2, 2> fullImage = ndarray::allocate(
        ndarray::makeVector(size, size)
    );
    ndarray::Array<double, 1, 1> compressedImage = ndarray::allocate(
        ndarray::makeVector(windowedFp.getNpix())
    );
    BOOST_CHECK_NO_THROW(windowedFp.compress(fullImage, compressedImage));

    ndarray::Array<double, 3, 3> matrix3d = ndarray::allocate(
        ndarray::makeVector(2, size, size)
    );
    ndarray::Array<double, 2, 2> matrix2d = ndarray::allocate(
        ndarray::makeVector(2, windowedFp.getNpix())
    );
    BOOST_CHECK_NO_THROW(windowedFp.compress(matrix3d, matrix2d));
}

BOOST_AUTO_TEST_CASE(WindowedFootprintExpand) {
    geom::PointI min = geom::PointI::makeXY(3,4);
    int size = 4;

    detection::Footprint fp(
        image::BBox(image::PointI(min.getX(), min.getY()), size, size)
    );
    image::BBox imageBBox = fp.getBBox();

    geom::BoxI geomBBox = geom::BoxI(min, geom::ExtentI::makeXY(size, size/2));
    multifit::WindowedFootprint windowedFp(fp, geomBBox);


    BOOST_WARN(windowedFp.getNpix()>0);
    ndarray::Array<double, 2, 2> fullImage = ndarray::allocate(
        ndarray::makeVector(size, size)
    );
    ndarray::Array<double, 1, 1> compressedImage = ndarray::allocate(
        ndarray::makeVector(windowedFp.getNpix())
    );

    BOOST_CHECK_NO_THROW(windowedFp.expand(compressedImage, fullImage));

    ndarray::Array<double, 3, 3> matrix3d = ndarray::allocate(
        ndarray::makeVector(2, size, size)
    );
    ndarray::Array<double, 2, 2> matrix2d = ndarray::allocate(
        ndarray::makeVector(2, windowedFp.getNpix())
    );
    BOOST_CHECK_NO_THROW(windowedFp.expand(matrix2d, matrix3d));
}
