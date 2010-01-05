// -*- LSST-C++ -*-
#include <iostream>
#include <cmath>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE CompressFunctor

#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"


#include "lsst/afw/geom/Box.h"
#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/matrices.h"
#include "lsst/meas/multifit/footprintUtils.h"

using namespace std;
namespace multifit = lsst::meas::multifit;


BOOST_AUTO_TEST_CASE(MatrixBasic) {
    lsst::afw::image::MaskedImage<float> testImg(5,5);
    *testImg.getImage() = 5;
    *testImg.getVariance() = 1;

    ndarray::Array<multifit::Pixel, 1, 1> img, var;
    ndarray::shallow(img) = ndarray::allocate<multifit::Allocator>(ndarray::makeVector(25));
    ndarray::shallow(var) = ndarray::allocate<multifit::Allocator>(ndarray::makeVector(25));

    lsst::afw::detection::Footprint fp(
        lsst::afw::image::BBox(lsst::afw::image::PointI(0,0), 5, 5)
    );


    multifit::CompressFunctor<float, unsigned short, float> func(testImg, img, var);
    func.apply(fp);
    
    ndarray::Array<multifit::Pixel, 1,1>::iterator iImg(img.begin()), iVar(var.begin());
    for(int i = 0; i < 25; ++i, ++iImg, ++iVar) {
        BOOST_CHECK_CLOSE(*iImg, 5.0, 0.00001);
        BOOST_CHECK_CLOSE(*iVar, 1.0, 0.00001);
    }
}
