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
#define BOOST_TEST_MODULE evaluator

#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"
#include "lsst/meas/multifit/Grid.h"
#include "lsst/meas/multifit/Evaluator.h"
#include "lsst/meas/multifit/ShapeletModelBasis.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/detection.h"

namespace geom = lsst::afw::geom;
namespace detection = lsst::afw::detection;
namespace image = lsst::afw::image;
namespace multifit = lsst::meas::multifit;


void checkEvaluator(
    multifit::Evaluator const & eval,
    multifit::Footprint const & fp,
    int nCoefficient, 
    int nParameters,
    multifit::Pixel dataValue,
    multifit::Pixel weightValue
) {
    BOOST_CHECK_EQUAL(eval.getGrid()->filterCount, 1);
    BOOST_CHECK_EQUAL(eval.getGrid()->coefficientCount, nCoefficient);
    BOOST_CHECK_EQUAL(eval.getGrid()->pixelCount, fp.getArea());
    BOOST_CHECK_EQUAL(eval.getGrid()->parameterCount, nParameters);

    BOOST_CHECK_EQUAL(eval.getGrid()->frames.size(), 1);
    multifit::grid::Frame frame = *(eval.getGrid()->frames.begin());
    BOOST_CHECK_EQUAL(frame.pixelCount, fp.getArea());
    BOOST_CHECK_EQUAL(frame.pixelOffset, 0);
    BOOST_CHECK_EQUAL(frame.frameIndex, 0);
    BOOST_CHECK_EQUAL(frame.id, 0);
    BOOST_CHECK_EQUAL(frame.data.getSize<0>(), fp.getArea());
    for(int i =0; i < fp.getArea(); ++i) {
        BOOST_CHECK_CLOSE(frame.data[i], dataValue, 0.00001);
    }

    BOOST_CHECK_EQUAL(frame.weights.getSize<0>(), fp.getArea());
    for(int i =0; i < fp.getArea(); ++i) {
        BOOST_CHECK_CLOSE(frame.weights[i], weightValue, 0.00001);
    }
   
    BOOST_CHECK_EQUAL(eval.getGrid()->objects.size(), 1);
    multifit::grid::Object object= *(eval.getGrid()->objects.begin());
    BOOST_CHECK_EQUAL(object.coefficientOffset, 0);
    BOOST_CHECK_EQUAL(object.coefficientCount, nCoefficient);
    BOOST_CHECK_EQUAL(eval.getGrid()->objects.begin()->sources.begin(), eval.getGrid()->sources.begin());
    BOOST_CHECK_EQUAL(eval.getGrid()->objects.begin()->sources.end(), eval.getGrid()->sources.end());
    BOOST_CHECK_EQUAL(object.id, 0);

    BOOST_CHECK_EQUAL(eval.getGrid()->sources.size(), 1);
    BOOST_CHECK_EQUAL(&eval.getGrid()->sources.begin()->object, eval.getGrid()->objects.begin());
    BOOST_CHECK_EQUAL(&eval.getGrid()->sources.begin()->frame, eval.getGrid()->frames.begin());

}
BOOST_AUTO_TEST_CASE(StarSourceConstruction) {
    geom::BoxI bbox = geom::BoxI(geom::PointI(25, 40), geom::ExtentI(10, 30));
    image::Exposure<float>::Ptr exp(new image::Exposure<float>(bbox));
    *exp->getMaskedImage().getImage() = 1.0;
    *exp->getMaskedImage().getVariance() = 4.0;
    detection::Psf::Ptr psf = detection::createPsf("DoubleGaussian", 19, 19, 2.0, 1.0);
    exp->setPsf(psf);
    detection::Footprint::Ptr fp(new detection::Footprint(bbox));
    geom::Point2D point(30.5, 55.8);
    
    multifit::Evaluator::Ptr eval = multifit::Evaluator::make(exp, fp, point);
    checkEvaluator(*eval, *fp, 1, 0, 1.0, 0.5);

    eval = multifit::Evaluator::make(exp, fp, point, false, true);
    checkEvaluator(*eval, *fp, 1, 2, 1.0, 0.5);
}

BOOST_AUTO_TEST_CASE(GalaxySourceConstruction) {
    geom::Point2D point(30.5, 55.8);
    geom::ellipses::Axes axes(10, 30, 0);
    geom::ellipses::Ellipse ellipse(axes, point); 
    geom::BoxI bbox = geom::BoxI(ellipse.computeEnvelope());
    image::Exposure<float>::Ptr exp(new image::Exposure<float>(bbox));
    *exp->getMaskedImage().getImage() = 1.0;
    *exp->getMaskedImage().getVariance() = 4.0;
    detection::Psf::Ptr psf = detection::createPsf("DoubleGaussian", 19, 19, 2.0, 1.0);
    exp->setPsf(psf);
    multifit::Footprint::Ptr fp(new multifit::Footprint(ellipse));
    multifit::ModelBasis::Ptr basis = multifit::ShapeletModelBasis::make(5);

    multifit::Evaluator::Ptr eval = multifit::Evaluator::make(exp, fp, basis, ellipse);
    checkEvaluator(*eval, *fp, basis->getSize(), 0, 1.0, 0.5);


    eval = multifit::Evaluator::make(exp, fp, basis, ellipse, true, true, false);
    checkEvaluator(*eval, *fp, basis->getSize(), 3, 1.0, 0.5);

    eval = multifit::Evaluator::make(exp, fp, basis, ellipse, true, true, true);
    checkEvaluator(*eval, *fp, basis->getSize(), 5, 1.0, 0.5);
}

BOOST_AUTO_TEST_CASE(DefinitionConstruction) {
    multifit::Definition definition;

    geom::BoxI box0 = geom::BoxI(geom::PointI(25, 40), geom::ExtentI(10, 30));
    detection::Psf::Ptr psf0 = detection::createPsf("DoubleGaussian", 19, 19, 2.0, 1.0);
    detection::Footprint::Ptr fp0(new detection::Footprint(box0));
    geom::Point2D point0(30.5, 55.8);

    lsst::ndarray::Array<multifit::Pixel, 1, 1> data0 = lsst::ndarray::allocate(
        lsst::ndarray::makeVector(fp0->getArea())
    );
    lsst::ndarray::Array<multifit::Pixel, 1, 1> weight0 = lsst::ndarray::allocate(
        lsst::ndarray::makeVector(fp0->getArea())
    );
    data0.deep()=1.0;
    weight0.deep()=0.5;

    int nFrame=0;
    //frame 0
    definition.frames.insert(
        multifit::definition::Frame(
            nFrame, 0, lsst::meas::multifit::Wcs::Ptr(), 
            psf0, fp0, data0, weight0
        )
    );
    nFrame++;
    //frame 1
    //
    lsst::ndarray::Array<multifit::Pixel, 1, 1> data1 = lsst::ndarray::allocate(
        lsst::ndarray::makeVector(fp0->getArea())
    );
    lsst::ndarray::Array<multifit::Pixel, 1, 1> weight1 = lsst::ndarray::allocate(
        lsst::ndarray::makeVector(fp0->getArea())
    );
    data1.deep()=2.0;
    weight1.deep()=1.0;

    definition.frames.insert(
        multifit::definition::Frame(
            nFrame, 1, lsst::meas::multifit::Wcs::Ptr(), 
            psf0, fp0, data1, weight1
        )
    );
    nFrame++;
    geom::Point2D point1(35.0, 57.5);
    geom::ellipses::Axes axes(10, 30, 0);
    geom::ellipses::Ellipse ellipse(axes, point1); 
    detection::Psf::Ptr psf1 = detection::createPsf("DoubleGaussian", 19, 19, 2.0, 1.0);
    multifit::Footprint::Ptr fp1(new multifit::Footprint(ellipse));
    lsst::ndarray::Array<multifit::Pixel, 1, 1> data2 = lsst::ndarray::allocate(
        lsst::ndarray::makeVector(fp1->getArea())
    );
    lsst::ndarray::Array<multifit::Pixel, 1, 1> weight2 = lsst::ndarray::allocate(
        lsst::ndarray::makeVector(fp1->getArea())
    );
    data2.deep() = 4.0;
    weight2.deep() = 2.0;
    
    //frame 2
    definition.frames.insert(
        multifit::definition::Frame(
            nFrame, 0, lsst::meas::multifit::Wcs::Ptr(), 
            psf1, fp1, data2, weight2
        )
    );
    nFrame++;

    multifit::ModelBasis::Ptr basis = multifit::ShapeletModelBasis::make(5);

    int nObject =0;
    int nParameters = 0;
    definition.objects.insert(
        multifit::definition::Object::makeStar(
            nObject++, 
            point0,
            false,
            true
        )
    );
    nParameters +=2;
    definition.objects.insert(
        multifit::definition::Object::makeGalaxy(
            nObject++, 
            basis,
            ellipse,
            true,
            true,
            false
        )
    );
    nParameters+=3;

    multifit::Evaluator::Ptr eval = multifit::Evaluator::make(definition);

    BOOST_CHECK_EQUAL(eval->getGrid()->filterCount, 2);
    BOOST_CHECK_EQUAL(eval->getGrid()->coefficientCount, 1+basis->getSize());
    BOOST_CHECK_EQUAL(eval->getGrid()->pixelCount, 2*fp0->getArea()+fp1->getArea());
    BOOST_CHECK_EQUAL(eval->getGrid()->parameterCount, nParameters);

    BOOST_CHECK_EQUAL(eval->getGrid()->frames.size(), nFrame);

    int pixelOffset =0, frameIndex =0;
    for(multifit::grid::Frame const * i = eval->getGrid()->frames.begin();
        i != eval->getGrid()->frames.end(); ++i
    ) {        
        multifit::definition::Frame const & j = definition.frames[i->id];
        BOOST_CHECK_EQUAL(i->id, j.id);
        BOOST_CHECK_EQUAL(i->frameIndex, frameIndex);
        BOOST_CHECK_EQUAL(i->wcs, j.wcs);
        BOOST_CHECK_EQUAL(i->psf, j.psf);
        BOOST_CHECK_EQUAL(i->footprint, j.footprint);
        BOOST_CHECK_EQUAL(i->pixelOffset, pixelOffset);
        BOOST_CHECK_EQUAL(i->pixelCount, j.footprint->getArea());


        BOOST_CHECK_EQUAL(i->weights.getSize<0>(), i->footprint->getArea());
        BOOST_CHECK_EQUAL(i->data.getSize<0>(), i->footprint->getArea());
        for(int pix=0; pix < i->footprint->getArea(); ++pix) {
            BOOST_CHECK_CLOSE(i->data[pix], j.data[pix], 0.00001);
            BOOST_CHECK_CLOSE(i->weights[pix], j.weights[pix], 0.00001);
        }
        pixelOffset += i->pixelCount;
        ++frameIndex;
    }

    int coefficientOffset=0;
    BOOST_CHECK_EQUAL(eval->getGrid()->objects.size(), nObject);
    for(multifit::grid::Object const * i = eval->getGrid()->objects.begin();
        i != eval->getGrid()->objects.end(); ++i
    ) {
        multifit::definition::Object const & j = definition.objects[i->id];
        BOOST_CHECK_EQUAL(i->id, j.id);
        BOOST_CHECK_EQUAL(i->basis, j.basis);
        BOOST_CHECK_EQUAL(i->isVariable, j.isVariable);
        BOOST_CHECK_EQUAL(i->radiusFactor, j.radiusFactor);
        BOOST_CHECK_EQUAL(i->coefficientCount, (j.basis) ? j.basis->getSize() : 1);
        BOOST_CHECK_EQUAL(i->coefficientOffset, coefficientOffset);

        for(multifit::grid::Object::SourceArray::const_iterator j = i->sources.begin();
            j != i->sources.end(); ++j
        ){
            BOOST_CHECK_EQUAL(&j->object, i);
        }
        coefficientOffset += i->coefficientCount;
    }

    BOOST_CHECK_EQUAL(eval->getGrid()->sources.size(), nObject*nFrame);
    for(multifit::grid::Source const * source = eval->getGrid()->sources.begin();
        source != eval->getGrid()->sources.end(); ++source
    ) {
        multifit::grid::Object const & object = multifit::grid::find(eval->getGrid()->objects, source->object.id);
        multifit::grid::Frame const & frame = multifit::grid::find(eval->getGrid()->frames, source->frame.id);
        
        BOOST_CHECK_EQUAL(&object, &source->object);
        BOOST_CHECK_EQUAL(&frame, &source->frame);
    }

}
