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
#define BOOST_TEST_MODULE FootprintUtils

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

BOOST_AUTO_TEST_CASE(clipAndMaskFootprintTest) {
    typedef lsst::afw::image::Mask<lsst::afw::image::MaskPixel> Mask;
    Mask::Ptr testMask = boost::make_shared<Mask>(8,8);
    testMask->setXY0(1, 1);
    *testMask = 0;
    detection::Footprint fp(
        lsst::afw::image::BBox(lsst::afw::image::PointI(0,0), 10, 10)
    );
    
    lsst::afw::image::BBox bbox(fp.getBBox());
    BOOST_CHECK_EQUAL(fp.getNpix(), 10*10);
    BOOST_CHECK_EQUAL(bbox.getX0(), 0);
    BOOST_CHECK_EQUAL(bbox.getX1(), 9);
    BOOST_CHECK_EQUAL(bbox.getY0(), 0);
    BOOST_CHECK_EQUAL(bbox.getY1(), 9);

    CONST_PTR(detection::Footprint) fixedFp(
        multifit::clipAndMaskFootprint<lsst::afw::image::MaskPixel>(fp, *testMask)
    );
    
    lsst::afw::geom::BoxI fixedBBox(lsst::afw::geom::convertToGeom(fixedFp->getBBox()));
    BOOST_CHECK_EQUAL(fixedFp->getNpix(), 8*8);
    BOOST_CHECK_EQUAL(fixedBBox.getMinX(), 1);
    BOOST_CHECK_EQUAL(fixedBBox.getMinY(), 1);
    BOOST_CHECK_EQUAL(fixedBBox.getMaxX(), 8);
    BOOST_CHECK_EQUAL(fixedBBox.getMaxY(), 8);


    (*testMask)(4,4) = 1;
    fixedFp = multifit::clipAndMaskFootprint<lsst::afw::image::MaskPixel>(fp, *testMask);
    fixedBBox = lsst::afw::geom::convertToGeom(fixedFp->getBBox());
    BOOST_CHECK_EQUAL(fixedFp->getNpix(), 8*8 - 1);
    BOOST_CHECK_EQUAL(fixedBBox.getMinX(), 1);
    BOOST_CHECK_EQUAL(fixedBBox.getMinY(), 1);
    BOOST_CHECK_EQUAL(fixedBBox.getMaxX(), 8);
    BOOST_CHECK_EQUAL(fixedBBox.getMaxY(), 8);
}

BOOST_AUTO_TEST_CASE(compressImageTest) {
    lsst::afw::image::MaskedImage<float> testImg(5,5);
    *testImg.getImage() = 5;
    *testImg.getVariance() = 1;
    *testImg.getMask() = 0;

    ndarray::Array<multifit::Pixel, 1, 1> img, var;


    lsst::afw::detection::Footprint fp(
        lsst::afw::image::BBox(lsst::afw::image::PointI(0,0), 5, 5)
    );

    ndarray::shallow(img) = ndarray::allocate<multifit::Allocator>(ndarray::makeVector(fp.getNpix()));
    ndarray::shallow(var) = ndarray::allocate<multifit::Allocator>(ndarray::makeVector(fp.getNpix()));
    multifit::compressImage<float, unsigned short, float>(fp, testImg, img, var);
    
    ndarray::Array<multifit::Pixel, 1,1>::iterator iImg(img.begin()), iVar(var.begin());
    for(int i = 0; i < fp.getNpix(); ++i, ++iImg, ++iVar) {
        BOOST_CHECK_CLOSE(*iImg, 5.0, 0.00001);
        BOOST_CHECK_CLOSE(*iVar, 1.0, 0.00001);
    }

    testImg.setXY0(1, 1);
    lsst::afw::detection::Footprint fp2(
        lsst::afw::image::BBox(lsst::afw::image::PointI(1,1), 2, 2)
    );

    ndarray::shallow(img) = ndarray::allocate<multifit::Allocator>(ndarray::makeVector(fp2.getNpix()));
    ndarray::shallow(var) = ndarray::allocate<multifit::Allocator>(ndarray::makeVector(fp2.getNpix()));
    multifit::compressImage(fp2, testImg, img, var);
    iImg = img.begin(), iVar = var.begin();
    for(int i = 0; i < fp2.getNpix(); ++i, ++iImg, ++iVar) {
        BOOST_CHECK_CLOSE(*iImg, 5.0, 0.00001);
        BOOST_CHECK_CLOSE(*iVar, 1.0, 0.00001);
    }
}
