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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE PointSourceModelProjection

#include "boost/make_shared.hpp"
#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"

#include "ndarray.hpp"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/Utils.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/geom/deprecated.h"

#include "lsst/meas/multifit/core.h"

#include "lsst/meas/multifit/footprintUtils.h"
#include "lsst/meas/multifit/ModelProjection.h"
#include "lsst/meas/multifit/ComponentModel.h"
#include "lsst/meas/multifit/components/PointSourceMorphology.h"
#include "lsst/meas/multifit/ModelFactory.h"

namespace afwDet = lsst::afw::detection;
namespace afwGeom = lsst::afw::geom;
namespace afwImg = lsst::afw::image;
namespace multifit = lsst::meas::multifit;

BOOST_AUTO_TEST_CASE(DerivativeTest) {
    afwGeom::PointD centroid = afwGeom::PointD::make(45,45);

    CONST_PTR(afwDet::Psf) psf = afwDet::createPsf("DoubleGaussian", 23, 23, 2.0);

    double flux = 5.45;

    afwGeom::AffineTransform transform;

    multifit::Model::Ptr psModel = 
        multifit::ModelFactory::createPointSourceModel(flux, centroid);

    CONST_PTR(afwDet::Footprint) fp = psModel->computeProjectionFootprint(
        psf, transform
    );
    PTR(multifit::ModelProjection) proj = psModel->makeProjection(
        psf, transform, fp
    );
    multifit::ParameterVector params(psModel->getNonlinearParameters());

    ndarray::Array<multifit::Pixel, 2, 2> analytic = ndarray::allocate<multifit::Allocator>(
        ndarray::makeVector(params.size(), fp->getNpix())
    );
    ndarray::Array<multifit::Pixel , 1, 1> plus = ndarray::allocate<multifit::Allocator>(
        ndarray::makeVector(fp->getNpix())
    );
    ndarray::Array<multifit::Pixel , 1, 1> minus = ndarray::allocate<multifit::Allocator>(
        ndarray::makeVector(fp->getNpix())
    );

    analytic = proj->computeNonlinearParameterDerivative();
    Eigen::MatrixXd numeric(fp->getNpix(), params.size());

    double eps = 2.2e-16;
    for(int i =0; i < params.size(); ++i) {
        double h = sqrt(eps)*params(i);
        params(i) += h;
        psModel->setNonlinearParameters(params);
        plus = proj->computeModelImage();

        params(i) -= 2*h;
        psModel->setNonlinearParameters(params);
        minus = proj->computeModelImage();

        Eigen::Map<Eigen::VectorXd> plusMap(plus.getData(), fp->getNpix(), 1);
        Eigen::Map<Eigen::VectorXd> minusMap(minus.getData(), fp->getNpix(), 1);
        numeric.col(i) = (plusMap-minusMap)/(2*h);
    }

    std::cerr << "analytic" << analytic;
    std::cerr << "numeric" << numeric;
}
