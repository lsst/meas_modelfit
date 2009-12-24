// -*- LSST-C++ -*-
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
    int nParam = 5;

    //ndarray is params by pixels
    ndarray::Array<Pixel, 2, 2> mArray = ndarray::allocate<multifit::Allocator>(
        ndarray::makeVector(nParam, nPix)
    );
    //but eigen is pix by params
    multifit::MatrixMapBlock mEigen = multifit::getMatrixView(mArray);

    //check casting from fully contiguous ndarray to eigen
    BOOST_CHECK_EQUAL(mEigen.cols(), nParam);
    BOOST_CHECK_EQUAL(mEigen.rows(), nPix);
    BOOST_CHECK_EQUAL(mEigen.cols(), mArray.getSize<0>());
    BOOST_CHECK_EQUAL(mEigen.rows(), mArray.getSize<1>());
    BOOST_CHECK_EQUAL(mEigen.data(), mArray.getData());



    int pixStart = 4, pixEnd = 9;
    int paramStart = 0, paramEnd = 4;
    //x== pix y == param
    lsst::afw::geom::BoxI mViewBox(
        lsst::afw::geom::PointI::makeXY(pixStart, paramStart),
        lsst::afw::geom::PointI::makeXY(pixEnd, paramEnd)
    );

    //now make a view into the master array
    ndarray::Array<Pixel, 2, 1> mView;
    ndarray::shallow(mView) = multifit::window(mArray, mViewBox);

    //and test casting from non fully contiguous ndarray to eigen
    multifit::MatrixMapBlock mViewEigen = multifit::getMatrixView(mView);
    BOOST_CHECK_EQUAL(mViewEigen.cols(), paramEnd - paramStart + 1);
    BOOST_CHECK_EQUAL(mViewEigen.rows(), pixEnd - pixStart + 1);
    BOOST_CHECK_EQUAL(mViewEigen.cols(), mView.getSize<0>());
    BOOST_CHECK_EQUAL(mViewEigen.rows(), mView.getSize<1>());
    BOOST_CHECK_EQUAL(mArray.getOwner(), mView.getOwner());
    BOOST_CHECK_EQUAL(mViewEigen.data(), mView.getData());



    int pixel = 8, param = 3;
    Pixel val = 42;
    mArray[param][pixel] = val;

    BOOST_CHECK_EQUAL(mEigen(pixel, param), val);
    BOOST_CHECK_EQUAL(mView[param - paramStart][pixel - pixStart], val);
    BOOST_CHECK_EQUAL(mViewEigen(pixel - pixStart, param -paramStart), val);
    BOOST_CHECK_NO_THROW(mView[paramEnd-paramStart][pixEnd-pixStart] = 0);
    BOOST_CHECK_NO_THROW(mViewEigen(pixEnd-pixStart, paramEnd-paramStart) = 0);
}
