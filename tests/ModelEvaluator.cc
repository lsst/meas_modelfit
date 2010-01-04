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
    typedef multifit::ModelEvaluator::Traits<float, lsst::afw::image::MaskPixel,
        lsst::afw::image::VariancePixel> Traits;

    multifit::ParameterVector linear(1), nonlinear (2);
    linear << 5;
    nonlinear << 25, 25;
    multifit::ModelFactory::ConstPtr factory = 
        multifit::ModelFactory::lookupFactory("PointSourceModelFactory");    
    multifit::Model::Ptr model = 
        factory->makeModel(1, linear.data(), nonlinear.data());

    multifit::ModelEvaluator evaluator(model);
    Traits::Exposure::Ptr exposure;

    multifit::PsfConstPtr psf = 
        measAlg::createPSF("DoubleGaussian", 19, 19, 1.5);

    Traits::CalibratedExposureList exposureList;

    //one exposure with full coverage
    exposure = boost::make_shared<Traits::Exposure>(50, 50);
    exposureList.push_back(Traits::CalibratedExposure(exposure, psf));
    
    //one exposure with partial coverage
    exposure = boost::make_shared<Traits::Exposure>(50, 25),
    exposureList.push_back(Traits::CalibratedExposure(exposure, psf));
    //one exposure with no coverage, by shifting image origin
    exposure = boost::make_shared<Traits::Exposure>(5, 5);
    exposure->getMaskedImage().setXY0(150, 150);
    exposureList.push_back(Traits::CalibratedExposure(exposure, psf));

    evaluator.setExposureList<float, lsst::afw::image::MaskPixel, 
        lsst::afw::image::VariancePixel>(exposureList);

    ndarray::Array<multifit::Pixel const, 1, 1> img;
    ndarray::Array<multifit::Pixel const, 1, 1> var;
    ndarray::Array<multifit::Pixel const, 1, 1> modelImg; 
    ndarray::Array<multifit::Pixel const, 2, 2> lpd, npd;

    ndarray::shallow(img) = evaluator.getImageVector();
    ndarray::shallow(var) = evaluator.getVarianceVector();
    ndarray::shallow(modelImg) = evaluator.computeModelImage();
    ndarray::shallow(lpd) = evaluator.computeLinearParameterDerivative();
    ndarray::shallow(npd) = evaluator.computeNonlinearParameterDerivative();
}
