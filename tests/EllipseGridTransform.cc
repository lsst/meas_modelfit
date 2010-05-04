
// -*- LSST-C++ -*-
#include <iostream>
#include <cmath>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE EllipseGridTransform

#include "boost/make_shared.hpp"
#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"

#include "ndarray.hpp"
#include "lsst/afw/image/Image.h"
#include "lsst/afw/image/Utils.h"
#include "lsst/afw/geom.h"
#include "lsst/meas/multifit/EllipseGridTransform.h"

using namespace std;
namespace multifit = lsst::meas::multifit;
namespace geom = lsst::afw::geom;

BOOST_AUTO_TEST_CASE(EllipseGridTransform) {
    lsst::afw::geom::ellipses::Axes axes(3, 1, 0);
    lsst::afw::geom::ExtentI dim = lsst::afw::geom::makeExtentI(120, 120);
    multifit::EllipseGridTransform egt(axes, dim);
    lsst::afw::geom::LinearTransform transform = egt;
    lsst::afw::image::Image<double> img(dim.getX(), dim.getY());
    for(int y = 0; y < img.getHeight(); ++y) {
        for( int x = 0; x < img.getWidth(); ++x) {
            lsst::afw::geom::Point2D p = lsst::afw::geom::makePointD(
                img.getX0() + x, 
                img.getY0() + ((y > img.getHeight()/2) ? (y-img.getHeight()) : y)
            );
            img(x, y) = transform(p).asVector().norm();
        }
    }

    img.writeFits("ellipseGridTest.fits");
}

