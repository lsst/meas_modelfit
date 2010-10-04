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
#define BOOST_TEST_MODULE ModelEvaluator

#include "boost/pointer_cast.hpp"
#include "boost/make_shared.hpp"
#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"

#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/coord/Coord.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/geom/deprecated.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/ComponentModel.h"
#include "lsst/meas/multifit/components/Astrometry.h"
#include "lsst/meas/multifit/components/PointSourceMorphology.h"
#include "lsst/meas/multifit/components/SersicMorphology.h"
#include "lsst/meas/multifit/ModelEvaluator.h"
#include "lsst/meas/multifit/ModelFactory.h"
#include "lsst/meas/multifit/RobustSersicCache.h"

namespace det = lsst::afw::detection;
namespace math = lsst::afw::math;
namespace image = lsst::afw::image;
namespace multifit = lsst::meas::multifit;
namespace geom = lsst::afw::geom;
namespace coord = lsst::afw::coord;

image::Wcs::Ptr makeWcs() {
    geom::PointD crPix= geom::makePointD(0,0);
    Eigen::Matrix2d cdMatrix;
    cdMatrix << 0.0001, 0.0, 0.0, 0.001;
    return boost::make_shared<image::Wcs>(geom::makePointD(45,45), crPix, cdMatrix);
}

image::MaskedImage<float> makeMaskedImage(
    multifit::Model::Ptr model,
    PTR(det::Psf) psf,
    geom::AffineTransform const & transform
) {

    det::Footprint::Ptr fp(model->computeProjectionFootprint(psf, transform));
    image::BBox bbox = fp->getBBox();


    std::cerr <<"computed footprint\n";

    image::MaskedImage<float> mi(
        bbox.getWidth(), 
        bbox.getHeight() 
    );    
    mi.setXY0(bbox.getX0(), bbox.getY0());


    multifit::ModelProjection::Ptr projection(model->makeProjection(psf, transform, fp));

    ndarray::Array<multifit::Pixel const, 1, 1> modelImage(projection->computeModelImage());
    ndarray::Array<multifit::Pixel, 1 ,1> variance(
        ndarray::allocate(ndarray::makeVector(fp->getNpix()))
    );

    variance = 0.25;
    multifit::expandImage(*fp, mi, modelImage, variance);

    return mi;
}

image::Exposure<float> makeExposure(
    multifit::Model::Ptr model, lsst::afw::geom::AffineTransform pixelToSky
) {
    PTR(det::Psf) psf = det::createPsf("DoubleGaussian", 9, 9, 1.5);
    PTR(image::Wcs) wcs = makeWcs();
    geom::AffineTransform transform = wcs->linearizeSkyToPixel(
        model->getPosition()
    );
    geom::AffineTransform pixToPix = pixelToSky*transform;
    std::cerr << "pix-to-pix: <" <<pixToPix[0] ;
    for(int i =1; i < 6; ++i) {
        std::cerr << ", " << pixToPix[i];        
    }
    std::cerr << ">\n";
    
    image::MaskedImage<float> mi = makeMaskedImage(model, psf, pixToPix);
    image::Exposure<float> exp(mi, *wcs);
    exp.setPsf(psf);

    return exp;
}

BOOST_AUTO_TEST_CASE(ConstructWithTransform) {
    double flux = 1;
    geom::Point2D pixel =geom::makePointD(45,45);

    geom::AffineTransform transform;
    multifit::Model::Ptr model = multifit::ModelFactory::createPointSourceModel(
        flux, 
        pixel
    );
    

    PTR(det::Psf) psf = det::createPsf("DoubleGaussian", 9,9, 1.5);

    multifit::ModelEvaluator eval(model, geom::AffineTransform());
    eval.setData<image::MaskedImage<float> >(makeMaskedImage(model, psf, transform), psf, transform);
}


BOOST_AUTO_TEST_CASE(PsModel) {
    double flux = 1;
    geom::Point2D centroid = geom::makePointD(45,45);

    geom::AffineTransform reference = geom::AffineTransform();
    multifit::Model::Ptr model = 
        multifit::ModelFactory::createPointSourceModel(
            flux, 
            centroid
        );

    std::list<image::Exposure<float> > exposureList;
    exposureList.push_back(makeExposure(model, reference));
    
    exposureList.push_back(makeExposure(model, reference));

    exposureList.push_back(makeExposure(model, reference));

    multifit::ModelEvaluator evaluator(model, reference);
    evaluator.setExposures(exposureList);

    BOOST_CHECK_EQUAL(evaluator.getNProjections(), 3);    
    BOOST_CHECK(evaluator.getNPixels() > 0);

    ndarray::Array<multifit::Pixel const, 1, 1> img;
    Eigen::Matrix<multifit::Pixel, Eigen::Dynamic, 1> modelImage;
    Eigen::Matrix<multifit::Pixel, Eigen::Dynamic, Eigen::Dynamic> lpd, npd;

    ndarray::shallow(img) = evaluator.getDataVector();
    modelImage = evaluator.computeModelImage();
    lpd = evaluator.computeLinearParameterDerivative();
    npd = evaluator.computeNonlinearParameterDerivative();
    
    //test for nan's in matrices
    for (int i = 0; i < evaluator.getNPixels(); ++i){
        BOOST_CHECK_EQUAL(img[i], img[i]);
        BOOST_CHECK_EQUAL(modelImage[i], modelImage[i]);

        for (int j = 0; j < evaluator.getLinearParameterSize(); ++j)
            BOOST_CHECK_EQUAL(lpd(i, j), lpd(i,j));

        for (int j = 0; j < evaluator.getNonlinearParameterSize(); ++j)
            BOOST_CHECK_EQUAL(npd(i,j), npd(i,j));
    }


}

BOOST_AUTO_TEST_CASE(SersicModel) {
    std::cerr << "SersicModel\n";
    //define the ellipse parameters in pixel coordinates
    double flux = 1;
    geom::Point2D pixel = geom::makePointD(45,45);
    geom::ellipses::Axes axes(25,30,0);

    multifit::Cache::ConstPtr cache;

    try {
        cache = multifit::Cache::load("testRobustCache", "RobustSersic", false);
    } catch (...){
        lsst::pex::policy::DefaultPolicyFile file("meas_multifit", "RobustSersicCache.paf", "tests");
        lsst::pex::policy::Policy pol;
        file.load(pol);
        cache = multifit::makeRobustSersicCache(pol);
    }
    multifit::components::SersicMorphology::setSersicCache(cache);

    //transform the ellipse parameters to be in sky coordinates
    geom::AffineTransform transform = geom::AffineTransform();
    std::cerr << "created cache\n";
    multifit::Model::Ptr model = multifit::ModelFactory::createSersicModel(
        flux, 
        pixel,
        axes, 
        1.0
    );
    std::cerr << "created model\n";
    
    multifit::ComponentModel::Ptr cm = boost::dynamic_pointer_cast<multifit::ComponentModel>(model);
    std::cerr << cm->getMorphology()->beginNonlinear() << "\n";
    std::cerr << model->getNonlinearParameterIter()+2 << "\n";
    double const * iter = model->getNonlinearParameterIter();
    for(int i =0; i < model->getNonlinearParameterSize(); ++i, ++iter){
        std::cerr << ", " << *iter;
    }
    std::cerr <<"\n";

    PTR(det::Psf) psf = det::createPsf("DoubleGaussian", 9,9, 1.5);
    image::MaskedImage<float> mi = makeMaskedImage(model, psf, transform);

    std::cerr << "created exposures\n";

    multifit::ModelEvaluator evaluator(model, transform);
    evaluator.setData(mi, psf, transform);


    BOOST_CHECK_EQUAL(evaluator.getNProjections(), 1);    
    BOOST_CHECK(evaluator.getNPixels() > 0);
    ndarray::Array<multifit::Pixel const, 1, 1> img;
    Eigen::Matrix<multifit::Pixel, Eigen::Dynamic, 1> modelImage;
    Eigen::Matrix<multifit::Pixel, Eigen::Dynamic, Eigen::Dynamic> lpd, npd;

    ndarray::shallow(img) = evaluator.getDataVector();
    modelImage = evaluator.computeModelImage();
    lpd = evaluator.computeLinearParameterDerivative();
    npd = evaluator.computeNonlinearParameterDerivative();
    
    //test for nan's in matrices
    for (int i = 0; i < evaluator.getNPixels(); ++i){
        BOOST_CHECK_EQUAL(img[i], img[i]);
        BOOST_CHECK_EQUAL(modelImage[i], modelImage[i]);

        for (int j = 0; j < evaluator.getLinearParameterSize(); ++j)
            BOOST_CHECK_EQUAL(lpd(i, j), lpd(i,j));

        for (int j = 0; j < evaluator.getNonlinearParameterSize(); ++j)
            BOOST_CHECK_EQUAL(npd(i,j), npd(i,j));
    }

    
}

