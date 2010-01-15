// -*- LSST-C++ -*-
#include <iostream>
#include <cmath>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE CompressFunctor

#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"


#include "lsst/afw/geom/Box.h"
#include "lsst/afw/image/Mask.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/image/Utils.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/matrices.h"
#include "lsst/meas/multifit/footprintUtils.h"

using namespace std;
namespace multifit = lsst::meas::multifit;
namespace detection = lsst::afw::detection;

BOOST_AUTO_TEST_CASE(clipAndMaskFootprint) {
    typedef lsst::afw::image::Mask<lsst::afw::image::MaskPixel> Mask;
    Mask::Ptr testMask = boost::make_shared<Mask>(8,8);
    testMask->setXY0(1, 1);
    *testMask = 0;
    multifit::FootprintConstPtr fp = boost::make_shared<detection::Footprint>(
        lsst::afw::image::BBox(lsst::afw::image::PointI(0,0), 10, 10)
    );
    
    lsst::afw::image::BBox bbox(fp->getBBox());
    BOOST_CHECK_EQUAL(fp->getNpix(), 10*10);
    BOOST_CHECK_EQUAL(bbox.getX0(), 0);
    BOOST_CHECK_EQUAL(bbox.getX1(), 9);
    BOOST_CHECK_EQUAL(bbox.getY0(), 0);
    BOOST_CHECK_EQUAL(bbox.getY1(), 9);

    multifit::FootprintConstPtr fixedFp(
        multifit::clipAndMaskFootprint<lsst::afw::image::MaskPixel>(fp, testMask)
    );
    
    lsst::afw::geom::BoxI fixedBBox(lsst::afw::geom::convertToGeom(fixedFp->getBBox()));
    BOOST_CHECK_EQUAL(fixedFp->getNpix(), 8*8);
    BOOST_CHECK_EQUAL(fixedBBox.getMinX(), 1);
    BOOST_CHECK_EQUAL(fixedBBox.getMinY(), 1);
    BOOST_CHECK_EQUAL(fixedBBox.getMaxX(), 8);
    BOOST_CHECK_EQUAL(fixedBBox.getMaxY(), 8);
}

BOOST_AUTO_TEST_CASE(CompressFunctorBasic) {
    lsst::afw::image::MaskedImage<float> testImg(5,5);
    *testImg.getImage() = 5;
    *testImg.getVariance() = 1;
    *testImg.getMask() = 0;

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

    testImg.setXY0(1, 1);
    lsst::afw::detection::Footprint fp2(
        lsst::afw::image::BBox(lsst::afw::image::PointI(1,1), 5, 5)
    );

    img = 0;
    var = 0;
    multifit::compressImage(fp2, testImg, img, var);
    iImg = img.begin(), iVar = var.begin();
    for(int i = 0; i < 25; ++i, ++iImg, ++iVar) {
        BOOST_CHECK_CLOSE(*iImg, 5.0, 0.00001);
        BOOST_CHECK_CLOSE(*iVar, 1.0, 0.00001);
    }
}
