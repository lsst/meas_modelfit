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
#define BOOST_TEST_MODULE matrices

#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"


#include "lsst/afw/geom/Box.h"
#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/matrices.h"

using namespace std;
namespace multifit = lsst::meas::multifit;

typedef multifit::Pixel Pixel;

BOOST_AUTO_TEST_CASE(MatrixBasic) {
    int nPix = 30;
    int nParam = 7;

    //ndarray is params by pixels
    ndarray::Array<Pixel, 2, 2> mArray = ndarray::allocate<multifit::Allocator>(
        ndarray::makeVector(nParam, nPix)
    );
    mArray = 1;
    //but eigen is pix by params
    multifit::MatrixMap mEigenMap(mArray.getData(), mArray.getStride<0>(), mArray.getSize<0>());
    multifit::MatrixMapBlock mEigen(mEigenMap, 0, 0, mArray.getSize<1>(), mArray.getSize<0>());

    //check casting from fully contiguous ndarray to eigen
    BOOST_CHECK_EQUAL(mEigen.cols(), nParam);
    BOOST_CHECK_EQUAL(mEigen.rows(), nPix);
    BOOST_CHECK_EQUAL(mEigen.cols(), mArray.getSize<0>());
    BOOST_CHECK_EQUAL(mEigen.rows(), mArray.getSize<1>());
    BOOST_CHECK_EQUAL(mEigen.data(), mArray.getData());
    for(int i = 0; i < mEigen.rows(); ++i) {
        for(int j = 0; j < mEigen.cols(); ++j) {
            BOOST_CHECK_EQUAL(mEigen(i,j), 1);
        }
    }
    int pixStart = 4, pixEnd = 9;
    int paramStart = 0, paramEnd = 7;

    //now make a view into the master array
    ndarray::Array<Pixel, 2, 1> mView(mArray[ndarray::view()(pixStart, pixEnd)]);
    mView = 2;

    //and test casting from non fully contiguous ndarray to eigen
    multifit::MatrixMap mViewEigenMap1(mView.getData(), mView.getStride<0>(), mView.getSize<0>());
    multifit::MatrixMapBlock mViewEigen1(mViewEigenMap1, 0, 0, mView.getSize<1>(), mView.getSize<0>());
    BOOST_CHECK_EQUAL(mViewEigen1.cols(), paramEnd - paramStart);
    BOOST_CHECK_EQUAL(mViewEigen1.rows(), pixEnd - pixStart);
    BOOST_CHECK_EQUAL(mViewEigen1.cols(), mView.getSize<0>());
    BOOST_CHECK_EQUAL(mViewEigen1.rows(), mView.getSize<1>());
    BOOST_CHECK_EQUAL(mArray.getOwner(), mView.getOwner());
    BOOST_CHECK_EQUAL(mViewEigen1.data(), mView.getData());
    for(int i = 0; i < mViewEigen1.rows(); ++i) {
        for(int j = 0; j < mViewEigen1.cols(); ++j) {
            BOOST_CHECK_EQUAL(mViewEigen1(i,j), 2);
        }
    }

    pixStart = 4;
    pixEnd = 9;
    paramStart = 1;
    paramEnd = 3;

    //now make a view into the master array
    ndarray::shallow(mView) = mArray[ndarray::view(paramStart, paramEnd)(pixStart, pixEnd)];
    mView = 3;

    //and test casting from non fully contiguous ndarray to eigen
    multifit::MatrixMap mViewEigenMap2(mView.getData(), mView.getStride<0>(), mView.getSize<0>());
    multifit::MatrixMapBlock mViewEigen2(mViewEigenMap2, 0, 0, mView.getSize<1>(), mView.getSize<0>());
    BOOST_CHECK_EQUAL(mViewEigen2.cols(), paramEnd - paramStart);
    BOOST_CHECK_EQUAL(mViewEigen2.rows(), pixEnd - pixStart);
    BOOST_CHECK_EQUAL(mViewEigen2.cols(), mView.getSize<0>());
    BOOST_CHECK_EQUAL(mViewEigen2.rows(), mView.getSize<1>());
    BOOST_CHECK_EQUAL(mArray.getOwner(), mView.getOwner());
    BOOST_CHECK_EQUAL(mViewEigen2.data(), mView.getData());
    for(int i = 0; i < mViewEigen2.rows(); ++i) {
        for(int j = 0; j < mViewEigen2.cols(); ++j) {
            BOOST_CHECK_EQUAL(mViewEigen2(i,j), 3);
        }
    }
    
    int pixel = 8, param = 2;
    Pixel val = 42;
    mArray[param][pixel] = val;

    BOOST_CHECK_EQUAL(&mEigen(0,0), &mArray[0][0]);
    BOOST_CHECK_EQUAL(&mEigen(pixel,param), &mArray[param][pixel]);
    BOOST_CHECK_EQUAL(mEigen(pixel, param), val);
    BOOST_CHECK_EQUAL(mView[param - paramStart][pixel - pixStart], val);
    BOOST_CHECK_EQUAL(mViewEigen2(pixel - pixStart, param -paramStart), val);
    BOOST_CHECK_NO_THROW(mView[paramEnd-paramStart-1][pixEnd-pixStart-1] = 0);
    BOOST_CHECK_NO_THROW(mViewEigen2(pixEnd-pixStart-1, paramEnd-paramStart-1) = 0);
}
