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
    typedef multifit::ModelEvaluator::Traits<double, lsst::afw::image::MaskPixel,
        lsst::afw::image::VariancePixel> Traits;

    multifit::ParameterVector linear(1), nonlinear (2);
    linear << 5;
    nonlinear << 25, 25;
    multifit::Model::Ptr model = multifit::makeModel("PointSource", linear, nonlinear);

    multifit::ModelEvaluator evaluator(model);

    lsst::afw::image::MaskedImage<double> testImg(50,50);
    //add a bogus variance
    *testImg.getVariance() = 0.1;

    image::Wcs wcs(
        image::PointD(1,1), image::PointD(1,1), Eigen::Matrix2d::Identity()
    );

    image::Exposure<double>::Ptr exposure = image::Exposure<double>::Ptr(
        new image::Exposure<double>(testImg, wcs)
    );
    multifit::PsfConstPtr psf = measAlg::createPSF("DoubleGaussian", 9, 9, 3);

    lsst::afw::image::Image<double> subImage(
        *testImg.getImage(), 
        lsst::afw::image::BBox(lsst::afw::image::PointI(20, 20), 9,9)
    );
    psf->getKernel()->computeImage(subImage, false);

    Traits::CalibratedExposureList exposureList;

    for(int i=0; i < 15; ++i) {
        exposureList.push_back(Traits::CalibratedExposure(exposure, psf));
    }
    evaluator.setExposureList<double, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>(exposureList);

    lsst::pex::policy::Policy::Ptr fitterPolicy(new lsst::pex::policy::Policy());
    fitterPolicy->add("terminationType", "iteration");    
    fitterPolicy->set("iterationMax", 2);
    
    multifit::SingleLinearParameterFitter fitter(fitterPolicy);
   
    multifit::SingleLinearParameterFitter::Result::Ptr result0 = fitter.apply(evaluator);

    multifit::SingleLinearParameterFitter::Result::Ptr result1 = fitter.apply(evaluator);
}
