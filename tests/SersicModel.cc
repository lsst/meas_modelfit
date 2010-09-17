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
#define BOOST_TEST_MODULE ModelProjection

#include "boost/make_shared.hpp"
#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"

#include "ndarray.hpp"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/Utils.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/geom/deprecated.h"

#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/Model.h"

#include "lsst/meas/multifit/footprintUtils.h"
#include "lsst/meas/multifit/ModelProjection.h"
#include "lsst/meas/multifit/ComponentModel.h"
#include "lsst/meas/multifit/FourierModelProjection.h"
#include "lsst/meas/multifit/components/SersicMorphology.h"
#include "lsst/meas/multifit/components/Astrometry.h"
#include "lsst/meas/multifit/ModelFactory.h"

using namespace std;
namespace multifit = lsst::meas::multifit;
namespace components = lsst::meas::multifit::components;
namespace geom = lsst::afw::geom;
namespace detection = lsst::afw::detection;

BOOST_AUTO_TEST_CASE(SersicModelProjection) {
    multifit::Cache::ConstPtr cache;
    try {
        cache = multifit::Cache::load("testCache", "Sersic", false);
    } catch (...){
        lsst::pex::policy::Policy pol;
        cache = multifit::makeSersicCache(pol);
    }
    geom::BoxD bounds = cache->getParameterBounds();
    cerr << bounds.getMinX() << " " << bounds.getMinY()<< endl; 
    cerr << bounds.getMaxX() << " " << bounds.getMaxY()<< endl; 
    multifit::components::SersicMorphology::setSersicCache(cache);

    lsst::afw::geom::PointD centroid = geom::PointD::make(0,0);

    //define ellipse in pixel coordinates
    lsst::afw::geom::ellipses::Axes axes(3, 1, 1.3);

    CONST_PTR(image::Wcs) wcs = boost::make_shared<image::Wcs>( 
        centroid, 
        geom::makePointD(0,0), 
        Eigen::Matrix2d::Identity()
    );
    //transform ellipse to sky coordinates
    geom::AffineTransform transform(wcs->linearizePixelToSky(centroid));
    axes.transform(transform).inPlace();

    lsst::afw::geom::ellipses::LogShear logShear(axes);
    double flux = 5.45;
    double sersicIndex = 2.0;

    multifit::Model::Ptr sgModel = multifit::ModelFactory::createSersicModel(
        flux, centroid, axes, sersicIndex
    );

    BOOST_CHECK_EQUAL(sgModel->getLinearParameterSize(), 1);
    BOOST_CHECK_EQUAL(sgModel->getNonlinearParameterSize(), 6);


    CONST_PTR(detection::Psf) psf = detection::createPsf("DoubleGaussian", 7, 7, 1.0);
    CONST_PTR(detection::Footprint) fp = sgModel->computeProjectionFootprint(psf, wcs);
    
    BOOST_CHECK(fp->getNpix() > 0);
    multifit::ModelProjection::Ptr projection = sgModel->makeProjection(psf, wcs, fp);
    BOOST_CHECK_EQUAL(projection->getModel(), sgModel);
   
    multifit::ParameterVector linear(sgModel->getLinearParameters());
    multifit::ParameterVector nonlinear(sgModel->getNonlinearParameters());

    BOOST_CHECK_EQUAL(linear[0], flux);
    BOOST_CHECK_EQUAL(nonlinear[0], centroid[0]);
    BOOST_CHECK_EQUAL(nonlinear[1], centroid[1]);
    BOOST_CHECK_EQUAL(nonlinear[2], logShear[0]);
    BOOST_CHECK_EQUAL(nonlinear[3], logShear[1]);
    BOOST_CHECK_EQUAL(nonlinear[4], logShear[2]);
    BOOST_CHECK_EQUAL(nonlinear[5], sersicIndex);

    ndarray::Array<multifit::Pixel const, 1, 1> modelImg; 
    ndarray::Array<multifit::Pixel const, 2, 1> lpd, npd;

    ndarray::shallow(modelImg) = projection->computeModelImage();
    ndarray::shallow(lpd) = projection->computeLinearParameterDerivative();
    ndarray::shallow(npd) = projection->computeNonlinearParameterDerivative();

    BOOST_CHECK_NO_THROW(projection->computeModelImage());
    BOOST_CHECK_NO_THROW(projection->computeLinearParameterDerivative());
    BOOST_CHECK_NO_THROW(projection->computeNonlinearParameterDerivative());
    
    //test for nan's in matrices
    for (int i = 0; i < modelImg.getSize<0>(); ++i){
        BOOST_CHECK_EQUAL(modelImg[i], modelImg[i]);

        for (int j = 0; j < sgModel->getLinearParameterSize(); ++j)
            BOOST_CHECK_EQUAL(lpd[j][i], lpd[j][i]);

        for (int j = 0; j < sgModel->getNonlinearParameterSize(); ++j)
            BOOST_CHECK_EQUAL(npd[j][i], npd[j][i]);
    }

    std::cerr << "npd" << npd <<std::endl;
    std::cerr << "lpd" << lpd <<std::endl;

    lsst::afw::image::BBox fpBbox = fp->getBBox();
    lsst::afw::image::Exposure<double> modelImage(
        fpBbox.getWidth(), fpBbox.getHeight(), *wcs
    );
    lsst::afw::image::MaskedImage<double> mi = modelImage.getMaskedImage();
    mi.setXY0(fpBbox.getLLC());

    multifit::expandImage(
        *fp, mi, projection->computeModelImage(),
        projection->computeLinearParameterDerivative()[0]
    );
    lsst::afw::detection::setMaskFromFootprint<lsst::afw::image::MaskPixel>(
        modelImage.getMaskedImage().getMask().get(), *fp, 1
    );
    modelImage.writeFits("sersicModelTest");
}
