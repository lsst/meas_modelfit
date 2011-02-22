// -*- LSST-C++ -*-
/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010, 2011 LSST Corporation.
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

#include "lsst/meas/multifit/shapelets/ShapeletModelBasis.h"
#include "lsst/meas/multifit/shapelets/ConvolvedShapeletModelBasis.h"
#include "lsst/ndarray/eigen.h"

namespace mf = lsst::meas::multifit;
namespace mfShapelets = lsst::meas::multifit::shapelets;
namespace afwShapelets = lsst::afw::math::shapelets;

void mfShapelets::ShapeletModelBasis::_evaluate(
    lsst::ndarray::Array<double, 2, 1> const & matrix,
    PTR(Footprint) const & footprint,
    lsst::afw::geom::Ellipse const & ellipse
) const {
    afw::geom::Ellipse fullEllipse(ellipse);
    fullEllipse.scale(_scale);
    afwShapelets::detail::HermiteEvaluator shapeletEvaluator(_order);
    afw::geom::AffineTransform transform = fullEllipse.getGridTransform();
    mf::Footprint::SpanList::const_iterator const spanEnd = footprint->getSpans().end();
    mf::Footprint::SpanList::const_iterator spanIter = footprint->getSpans().begin();
    lsst::ndarray::Array<double, 2, 1>::Iterator pixelIter = matrix.begin();
    for (; spanIter != spanEnd; ++spanIter) {
        afw::detection::Span const & span = **spanIter;
        for (int x = span.getX0(); x <= span.getX1(); ++x) {
            shapeletEvaluator.fillEvaluation(*pixelIter, transform(afw::geom::Point2D(x, span.getY())));
        }
    }
}

mf::ModelBasis::Ptr mfShapelets::ShapeletModelBasis::_convolve(
    lsst::meas::multifit::LocalPsf::Ptr const & psf
) const {
    LocalPsf::Shapelet shapeletPsf = psf->asShapelet(afwShapelets::HERMITE);
    return boost::make_shared<ConvolvedShapeletModelBasis>(*this, shapeletPsf);
}
