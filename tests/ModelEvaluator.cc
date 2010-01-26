// -*- LSST-C++ -*-
#include <iostream>
#include <cmath>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ModelEvaluator

#include "boost/pointer_cast.hpp"
#include "boost/make_shared.hpp"
#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"

#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/image/MaskedImage.h"

#include "lsst/afw/geom/Box.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/PointSourceModelFactory.h"
#include "lsst/meas/multifit/ModelEvaluator.h"

namespace math = lsst::afw::math;
namespace image = lsst::afw::image;
namespace multifit = lsst::meas::multifit;
namespace geom = lsst::afw::geom;
namespace measAlg = lsst::meas::algorithms;

BOOST_AUTO_TEST_CASE(ModelBasic) {
    multifit::PointSourceModelFactory psFactory;
    multifit::Model::Ptr psModel = psFactory.makeModel(1, geom::makePointD(45, 45));

    image::PointD crPix(0, 0), crVal(45,45);
    Eigen::Matrix2d cdMatrix(Eigen::Matrix2d::Identity()*0.0001);
    multifit::Wcs::Ptr wcs = boost::make_shared<multifit::Wcs> (crVal, crPix, cdMatrix);

    multifit::Psf::Ptr psf = measAlg::createPSF("DoubleGaussian", 19, 19, 2);
    multifit::FootprintConstPtr fp(psModel->computeProjectionFootprint(psf, wcs));
    image::BBox bbox = fp->getBBox();

    image::MaskedImage<double> mi(
        bbox.getWidth(), 
        bbox.getHeight() 
    );
    mi.setXY0(bbox.getX0(), bbox.getY0());
    *(mi.getMask()) = 0;

    multifit::ModelProjection::Ptr projection(psModel->makeProjection(psf, wcs, fp));
    ndarray::Array<multifit::Pixel const, 1, 1> modelImage(projection->computeModelImage());
    ndarray::Array<multifit::Pixel, 1 ,1> variance(ndarray::allocate(ndarray::makeVector(fp->getNpix())));
    variance = 0.5*0.5;

    multifit::expandImage(*fp, mi, modelImage, variance);

    std::list<multifit::CharacterizedExposure<multifit::Pixel>::Ptr> exposureList;

    multifit::CharacterizedExposure<double>::Ptr toAdd;
    //one exposure with full coverage
    toAdd.reset(new multifit::CharacterizedExposure<multifit::Pixel>(mi, *wcs, psf));
    exposureList.push_back(toAdd);
    
    //one exposure with partial coverage
    image::BBox partial(image::PointI(0,0), bbox.getWidth() / 2, bbox.getHeight() /2);
    image::MaskedImage<double> sub(mi, partial);
    toAdd.reset(new multifit::CharacterizedExposure<multifit::Pixel> (sub, *wcs, psf));
    exposureList.push_back(toAdd);

    //one exposure with no coverage, by shifting image origin
    bbox.shift(-bbox.getX0(), -bbox.getY0());
    image::MaskedImage<double> offset(mi, bbox, true);
    offset.setXY0(-1000, -1000);
    toAdd.reset(new multifit::CharacterizedExposure<multifit::Pixel> (offset, *wcs, psf));
    exposureList.push_back(toAdd);

    multifit::ModelEvaluator evaluator(psModel, exposureList);

    BOOST_CHECK_EQUAL(evaluator.getNProjections(), 2);    
    BOOST_CHECK(evaluator.getNPixels() < fp->getNpix()*2);
    ndarray::Array<multifit::Pixel const, 1, 1> img;
    ndarray::Array<multifit::Pixel const, 1, 1> var;
    ndarray::Array<multifit::Pixel const, 1, 1> modelImg; 
    ndarray::Array<multifit::Pixel const, 2, 2> lpd, npd;

    ndarray::shallow(img) = evaluator.getDataVector();
    ndarray::shallow(var) = evaluator.getVarianceVector();
    ndarray::shallow(modelImg) = evaluator.computeModelImage();
    ndarray::shallow(lpd) = evaluator.computeLinearParameterDerivative();
    ndarray::shallow(npd) = evaluator.computeNonlinearParameterDerivative();
    
    for (int i = 0; i < evaluator.getNPixels(); ++i){
        BOOST_CHECK_EQUAL(img[i], img[i]);
        BOOST_CHECK_EQUAL(var[i], var[i]);
        BOOST_CHECK_EQUAL(modelImg[i], modelImg[i]);

        for (int j = 0; j < evaluator.getLinearParameterSize(); ++j)
            BOOST_CHECK_EQUAL(lpd[j][i], lpd[j][i]);

        for (int j = 0; j < evaluator.getNonlinearParameterSize(); ++j)
            BOOST_CHECK_EQUAL(npd[j][i], npd[j][i]);
    }
}
