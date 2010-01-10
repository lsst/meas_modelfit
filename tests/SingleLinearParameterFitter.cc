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

    multifit::ParameterVector linear(1), nonlinear (2);
    linear << 1;
    nonlinear << 25, 25;
    multifit::Model::Ptr model = multifit::makeModel("PointSource", linear, nonlinear);

    multifit::ModelEvaluator evaluator(model);

    image::Wcs wcs(
        image::PointD(0,0), image::PointD(0,0), Eigen::Matrix2d::Identity()
    );

    multifit::Psf::Ptr psf = measAlg::createPSF("DoubleGaussian", 9, 9, 3);
    multifit::CharacterizedExposure<double>::Ptr exposure = 
        boost::make_shared<multifit::CharacterizedExposure<double> >(50, 50, wcs, psf);

    //add a bogus variance
    lsst::afw::image::Image<float>::Ptr variance = exposure->getMaskedImage().getVariance();
    *variance = 0.1;

    lsst::afw::image::Image<double> subImage(
        *exposure->getMaskedImage().getImage(), 
        lsst::afw::image::BBox(lsst::afw::image::PointI(20, 20), 9,9)
    );
    psf->getKernel()->computeImage(subImage, true);

    std::list<multifit::CharacterizedExposure<double>::Ptr> exposureList;

    for(int i=0; i < 15; ++i) {
        exposureList.push_back(exposure);
    }
    evaluator.setExposureList(exposureList);

    lsst::pex::policy::Policy::Ptr fitterPolicy(new lsst::pex::policy::Policy());
    fitterPolicy->add("terminationType", "iteration");    
    fitterPolicy->set("iterationMax", 2);
    
    multifit::SingleLinearParameterFitter fitter(fitterPolicy);
   
    multifit::SingleLinearParameterFitter::Result::Ptr result0 = fitter.apply(evaluator);

    multifit::SingleLinearParameterFitter::Result::Ptr result1 = fitter.apply(evaluator);
}
