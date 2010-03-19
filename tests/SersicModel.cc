// -*- LSST-C++ -*-
#include <iostream>
#include <cmath>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ModelProjection

#include "boost/make_shared.hpp"
#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"

#include "ndarray.hpp"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/Utils.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/meas/algorithms/PSF.h"
#include "lsst/afw/geom/deprecated.h"

#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/Model.h"

#include "lsst/meas/multifit/footprintUtils.h"
#include "lsst/meas/multifit/ModelProjection.h"
#include "lsst/meas/multifit/FourierModelProjection.h"
#include "lsst/meas/multifit/components/SersicMorphology.h"
#include "lsst/meas/multifit/SersicModelFactory.h"

using namespace std;
namespace measAlg = lsst::meas::algorithms;
namespace multifit = lsst::meas::multifit;
namespace geom = lsst::afw::geom;
namespace detection = lsst::afw::detection;

BOOST_AUTO_TEST_CASE(SersicModelProjection) {
    multifit::SersicModelFactory sgFactory;
    lsst::afw::geom::PointD centroid = geom::PointD::make(0,0);
    lsst::afw::geom::ellipses::Axes::Ptr axes(
        new lsst::afw::geom::ellipses::Axes(3, 1, 1.3)
    );
    lsst::afw::geom::ellipses::LogShear logShear(*axes);
    double flux = 5.45;
    double sersicIndex = 2.0;

    multifit::Model::Ptr sgModel = sgFactory.makeModel(
        flux, 
        centroid, 
        axes, 
        sersicIndex
    );
    
    BOOST_CHECK_EQUAL(sgModel->getLinearParameterSize(), 1);
    BOOST_CHECK_EQUAL(sgModel->getNonlinearParameterSize(), 6);
    multifit::WcsConstPtr wcs = boost::make_shared<multifit::Wcs>( 
        lsst::afw::geom::convertToImage(centroid), 
        image::PointD(0,0), 
        Eigen::Matrix2d::Identity()
    );

    multifit::PsfConstPtr psf = measAlg::createPSF("DoubleGaussian", 7, 7, 1.0);
    multifit::FootprintConstPtr fp = sgModel->computeProjectionFootprint(psf, wcs);
    
    BOOST_CHECK(fp->getNpix() > 0);
    multifit::ModelProjection::Ptr projection = sgModel->makeProjection(psf, wcs, fp);
    BOOST_CHECK_EQUAL(projection->getModel(), sgModel);
   
    multifit::ParameterVector linear(sgModel->getLinearParameterSize());
    multifit::ParameterVector nonlinear(sgModel->getNonlinearParameterSize());

    linear << sgModel->getLinearParameters();
    nonlinear << sgModel->getNonlinearParameters();

    BOOST_CHECK_EQUAL(linear[0], flux);
    BOOST_CHECK_EQUAL(nonlinear[0], centroid[0]);
    BOOST_CHECK_EQUAL(nonlinear[1], centroid[1]);
    BOOST_CHECK_EQUAL(nonlinear[2], logShear[0]);
    BOOST_CHECK_EQUAL(nonlinear[3], logShear[1]);
    BOOST_CHECK_EQUAL(nonlinear[4], logShear[2]);
    BOOST_CHECK_EQUAL(nonlinear[5], sersicIndex);

    BOOST_CHECK_NO_THROW(projection->computeModelImage());
    BOOST_CHECK_NO_THROW(projection->computeLinearParameterDerivative());
    BOOST_CHECK_NO_THROW(projection->computeNonlinearParameterDerivative());

    lsst::afw::image::BBox fpBbox = fp->getBBox();
    lsst::afw::image::Exposure<double> modelImage(
        fpBbox.getWidth(), fpBbox.getHeight(), *wcs
    );
    modelImage.getMaskedImage().setXY0(fpBbox.getLLC());
    multifit::expandImage(
        *fp, modelImage.getMaskedImage(), projection->computeModelImage(),
        projection->computeLinearParameterDerivative()[0]
    );
    lsst::afw::detection::setMaskFromFootprint<lsst::afw::image::MaskPixel>(
        modelImage.getMaskedImage().getMask().get(), *fp, 1
    );
    modelImage.writeFits("sersicModelTest");
}
