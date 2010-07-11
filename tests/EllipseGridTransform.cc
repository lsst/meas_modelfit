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

