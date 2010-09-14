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
#define BOOST_TEST_MODULE MinuitFitter

#include "boost/pointer_cast.hpp"
#include "boost/make_shared.hpp"
#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"

#include "lsst/afw/image/Exposure.h"
#include <Eigen/Core>
#include <Eigen/Array>
#include "lsst/afw/math/Random.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/geom/deprecated.h"
#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/ModelEvaluator.h"
#include "lsst/meas/multifit/ComponentModel.h"
#include "lsst/meas/multifit/MinuitFitter.h"
#include "lsst/meas/multifit/components/Astrometry.h"
#include "lsst/meas/multifit/components/PointSourceMorphology.h"
#include "lsst/meas/multifit/ModelFactory.h"

namespace math = lsst::afw::math;
namespace image = lsst::afw::image;
namespace multifit = lsst::meas::multifit;
namespace geom = lsst::afw::geom;
namespace detection = lsst::afw::detection;

BOOST_AUTO_TEST_CASE(BasicFitter) {
    geom::PointD centroid = geom::PointD::make(35,65);
    double flux = 34.45;
    multifit::components::PointSourceMorphology::Ptr morphology =
        multifit::components::PointSourceMorphology::create(flux);
    multifit::components::Astrometry astrometry(centroid);

    //define ellipse in pixel coordinates
    geom::ellipses::Axes axes(30, 15, 1.3);

    geom::PointD crPix(0), crVal(centroid);
    Eigen::Matrix2d cdMatrix(Eigen::Matrix2d::Identity());
    image::Wcs::Ptr wcs = boost::make_shared<image::Wcs> (crVal, crPix, cdMatrix);

    //transform ellipse to sky coordinates
    geom::AffineTransform transform(wcs->linearizePixelToSky(centroid));
    axes.transform(transform).inPlace();
    
    multifit::Model::Ptr model = 
        multifit::ModelFactory::createExponentialModel(flux, centroid, axes);
    //    multifit::ModelFactory::createPointSourceModel(flux, centroid);
    detection::Psf::Ptr psf = detection::createPsf("DoubleGaussian", 19, 19, 2);

    CONST_PTR(detection::Footprint) fp(model->computeProjectionFootprint(psf, wcs));
    image::BBox bbox = fp->getBBox();
    
    image::Exposure<double>::Ptr exposure = 
        boost::make_shared<image::Exposure<double> >(bbox.getWidth(), bbox.getHeight(), *wcs);
    exposure->setPsf(psf);
    lsst::afw::image::MaskedImage<double> mi = exposure->getMaskedImage();
    mi.setXY0(bbox.getX0(), bbox.getY0());
    *mi.getMask() = 0;

    multifit::ModelProjection::Ptr projection(model->makeProjection(psf, wcs, fp));
    ndarray::Array<multifit::Pixel const, 1, 1> modelImage(projection->computeModelImage());
    ndarray::Array<multifit::Pixel, 1 ,1> variance(ndarray::allocate(ndarray::makeVector(fp->getNpix())));
    variance = 0.5*0.5;

    multifit::expandImage(*fp, mi, modelImage, variance);
   
    math::Random random;
    image::Image<double> randomImg(mi.getWidth(), mi.getHeight());
    math::randomUniformImage<image::Image<double> >(&randomImg, random);
    *mi.getImage() += randomImg;

    std::list<image::Exposure<double>::Ptr> exposureList;
    for(int i=0; i < 5; ++i) {
        exposureList.push_back(exposure);
    }
    multifit::ModelEvaluator::Ptr evaluator(new multifit::ModelEvaluator(model));
    evaluator->setExposureList<double, image::MaskPixel, image::VariancePixel>(
        exposureList
    );
       
    std::vector<double> errors(
        evaluator->getLinearParameterSize() + evaluator->getNonlinearParameterSize(),
        0.001
    );
    lsst::pex::policy::Policy::Ptr fitterPolicy(new lsst::pex::policy::Policy());
    std::cerr << "nPix: " << evaluator->getNPixels() << "\n";
    multifit::MinuitAnalyticFitter analytic(fitterPolicy);
    multifit::MinuitNumericFitter numeric(fitterPolicy);

    multifit::MinuitFitterResult numericResult = numeric.apply(evaluator, errors);
    multifit::MinuitFitterResult analyticResult = analytic.apply(evaluator, errors);

}


