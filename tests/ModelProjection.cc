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
#include "lsst/meas/algorithms/PSF.h"

#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/Model.h"

#include "lsst/meas/multifit/ModelProjection.h"
#include "lsst/meas/multifit/FourierModelProjection.h"
#include "lsst/meas/multifit/PointSourceModelFactory.h"

#include "lsst/meas/multifit/components/Astrometry.h"
#include "lsst/meas/multifit/components/PointSourceMorphology.h"

using namespace std;
namespace measAlg = lsst::meas::algorithms;
namespace multifit = lsst::meas::multifit;
namespace geom = lsst::afw::geom;
namespace detection = lsst::afw::detection;

BOOST_AUTO_TEST_CASE(FourierModelProjection) {
    multifit::PointSourceModelFactory psFactory;
    lsst::afw::geom::PointD centroid = geom::PointD::makeXY(35,65);
    double flux = 34.45;
    multifit::Model::Ptr psModel = psFactory.makeModel(flux, centroid);
    
    BOOST_CHECK_EQUAL(psModel->getLinearParameterSize(), 1);
    BOOST_CHECK_EQUAL(psModel->getNonlinearParameterSize(), 2);

    multifit::WcsConstPtr wcs = boost::make_shared<multifit::Wcs>();
    multifit::PsfConstPtr psf = measAlg::createPSF("DoubleGaussian", 7, 7, 1.5);
    multifit::FootprintConstPtr fp = psModel->computeProjectionFootprint(psf, wcs);
    
    BOOST_CHECK(fp->getNpix() > 0);
    multifit::ModelProjection::Ptr projection = psModel->makeProjection(psf, wcs, fp);
    BOOST_CHECK_EQUAL(projection->getModel(), psModel);
   
    multifit::ParameterVector linear(psModel->getLinearParameterSize());
    multifit::ParameterVector nonlinear(psModel->getNonlinearParameterSize());

    linear = psModel->getLinearParameterVector();
    nonlinear = psModel->getNonlinearParameterVector();

    BOOST_CHECK_EQUAL(linear[0], flux);
    BOOST_CHECK_EQUAL(nonlinear[0], centroid[0]);
    BOOST_CHECK_EQUAL(nonlinear[1], centroid[1]);

    multifit::FourierModelProjection::Ptr asFourierProjection =
        boost::static_pointer_cast<multifit::FourierModelProjection>(projection);

    BOOST_CHECK(asFourierProjection);
    BOOST_CHECK_NO_THROW(projection->computeModelImage());
    BOOST_CHECK_NO_THROW(projection->computeLinearParameterDerivative());
    BOOST_CHECK_NO_THROW(projection->computeNonlinearParameterDerivative());
}
