// -*- LSST-C++ -*-
#include <iostream>
#include <cmath>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SingleLinearParameterFitter

#include "boost/pointer_cast.hpp"
#include "boost/make_shared.hpp"
#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"

#include "lsst/afw/image/Exposure.h"
#include <Eigen/Core>
#include <Eigen/Array>
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/PointSourceModelFactory.h"
#include "lsst/meas/multifit/ModelEvaluator.h"
#include "lsst/meas/multifit/SingleLinearParameterFitter.h"

namespace math = lsst::afw::math;
namespace image = lsst::afw::image;
namespace multifit = lsst::meas::multifit;
namespace geom = lsst::afw::geom;
namespace measAlg = lsst::meas::algorithms;

BOOST_AUTO_TEST_CASE(FitterBasic) {
    multifit::PointSourceModelFactory psFactory;
    multifit::Model::Ptr psModel = psFactory.makeModel(1, geom::makePointD(45, 45));

    image::PointD crPix(0, 0), crVal(45,45);
    Eigen::Matrix2d cdMatrix(Eigen::Matrix2d::Identity()*0.0001);
    image::Wcs::Ptr wcs = boost::make_shared<image::Wcs> (crVal, crPix, cdMatrix);

    multifit::Psf::Ptr psf = measAlg::createPSF("DoubleGaussian", 19, 19, 2);
    multifit::FootprintConstPtr fp(psModel->computeProjectionFootprint(psf, wcs));
    image::BBox bbox = fp->getBBox();

    multifit::CharacterizedExposure<double>::Ptr exposure = 
        boost::make_shared<multifit::CharacterizedExposure<double> >(bbox.getWidth(), bbox.getHeight(), *wcs, psf);
    exposure->getMaskedImage().setXY0(bbox.getX0(), bbox.getY0());
    *exposure->getMaskedImage().getMask() = 0;

    multifit::ModelProjection::Ptr projection(psModel->makeProjection(psf, wcs, fp));
    ndarray::Array<multifit::Pixel const, 1, 1> modelImage(projection->computeModelImage());
    ndarray::Array<multifit::Pixel, 1 ,1> variance(ndarray::allocate(ndarray::makeVector(fp->getNpix())));
    variance = 0.5*0.5;

    multifit::expandImage(*fp, exposure->getMaskedImage(), modelImage, variance);

    std::list<multifit::CharacterizedExposure<double>::Ptr> exposureList;
    for(int i=0; i < 5; ++i) {
        exposureList.push_back(exposure);
    }
    multifit::ModelEvaluator evaluator(psModel, exposureList);
       
    lsst::pex::policy::Policy::Ptr fitterPolicy(new lsst::pex::policy::Policy());
    fitterPolicy->add("terminationType", "iteration");    
    fitterPolicy->add("terminationType", "dChisq");
    fitterPolicy->set("iterationMax", 5);
    fitterPolicy->set("dChisqThreshold", 0.001);
    
    multifit::SingleLinearParameterFitter fitter(fitterPolicy);
   
    multifit::SingleLinearParameterFitter::Result::Ptr result0 = fitter.apply(evaluator);

    multifit::SingleLinearParameterFitter::Result::Ptr result1 = fitter.apply(evaluator);
}
