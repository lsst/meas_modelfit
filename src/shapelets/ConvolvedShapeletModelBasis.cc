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

mfShapelets::ConvolvedShapeletModelBasis::ConvolvedShapeletModelBasis(
    ShapeletModelBasis const & basis,
    lsst::afw::math::shapelets::ShapeletFunction const & psf
) : ModelBasis(basis.getSize()),
    _convolution(boost::make_shared<afwShapelets::detail::HermiteConvolution>(basis.getOrder(), psf)),
    _frontBasis(ShapeletModelBasis::make(_convolution->getRowOrder(), 1.0)),
    _scale(basis.getScale())
{}

void mfShapelets::ConvolvedShapeletModelBasis::_evaluate(
    lsst::ndarray::Array<double,2,1> const & matrix,
    lsst::afw::detection::Footprint::Ptr const & footprint,
    lsst::afw::geom::Ellipse const & ellipse
) const {
    ndarray::Array<double,2,2> frontMatrix(
        ndarray::allocate(footprint->getArea(), _frontBasis->getSize())
    );
    afw::geom::Ellipse frontEllipse(ellipse);
    frontEllipse.scale(_scale);
    ndarray::Array<afwShapelets::Pixel const,2,2> convolutionMatrix = 
        _convolution->evaluate(frontEllipse.getCore());
    _frontBasis->evaluate(frontMatrix, footprint, frontEllipse);
    ndarray::viewAsEigen(matrix) = 
        ndarray::viewAsEigen(frontMatrix) * ndarray::viewAsEigen(convolutionMatrix) / _scale;
}
