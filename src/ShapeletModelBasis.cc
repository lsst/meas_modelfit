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

#include "lsst/meas/multifit/ShapeletModelBasis.h"
#include "lsst/meas/multifit/CompoundShapeletModelBasis.h"
#include "lsst/afw/math/shapelets/HermiteConvolution.h"
#include "lsst/ndarray/eigen.h"

namespace afwShapelets = lsst::afw::math::shapelets;

namespace lsst { namespace meas { namespace multifit {

ConvolvedShapeletModelBasis::ConvolvedShapeletModelBasis(
    ShapeletModelBasis const & basis,
    afwShapelets::ShapeletFunction const & psf
) : ModelBasis(basis.getSize()),
    _convolution(boost::make_shared<afwShapelets::HermiteConvolution>(basis.getOrder(), psf)),
    _frontBasis(ShapeletModelBasis::make(_convolution->getRowOrder(), 1.0)),
    _scale(basis.getScale())
{}

void ConvolvedShapeletModelBasis::_integrate(lsst::ndarray::Array<Pixel, 1, 1> const & vector) const {
    throw LSST_EXCEPT(
        lsst::pex::exceptions::LogicErrorException,
        "Cannot integrate convolved basis."
    );      
}

void ConvolvedShapeletModelBasis::_evaluate(
    ndarray::Array<Pixel, 2, 1> const & matrix,
    CONST_PTR(Footprint) const & footprint,
    afw::geom::Ellipse const & ellipse
) const {
    ndarray::Array<Pixel,2,2> frontMatrix(
        ndarray::allocate(footprint->getArea(), _frontBasis->getSize())
    );
    afw::geom::Ellipse frontEllipse(ellipse);
    frontEllipse.scale(_scale);
    ndarray::Array<afwShapelets::Pixel const,2,2> convolutionMatrix = 
        _convolution->evaluate(frontEllipse);
    _frontBasis->evaluate(frontMatrix, footprint, frontEllipse);
    ndarray::viewAsEigen(matrix) = 
        ndarray::viewAsEigen(frontMatrix) * ndarray::viewAsEigen(convolutionMatrix) / _scale;
}

Eigen::MatrixXd ConvolvedShapeletModelBasis::computeInnerProductMatrix(
    lsst::afw::geom::ellipses::BaseCore const & ellipse
) const {
    afw::geom::Ellipse frontEllipse(ellipse);
    frontEllipse.scale(_scale);
    ndarray::Array<Pixel const,2,2> c = _convolution->evaluate(frontEllipse);
    double r = frontEllipse.getCore().getDeterminantRadius();
    Eigen::MatrixXd m = afw::math::shapelets::HermiteEvaluator::computeInnerProductMatrix(
        _frontBasis->getOrder(), _frontBasis->getOrder(), 1.0 / r, 1.0 / r
    );
    m /= _scale * _scale;
    return ndarray::viewAsTransposedEigen(c) * m * ndarray::viewAsEigen(c);
}

int & ShapeletModelBasis::getPsfShapeletOrderRef() {
    static int v = 4;
    return v;
}

Eigen::MatrixXd ShapeletModelBasis::computeInnerProductMatrix(
    lsst::afw::geom::ellipses::BaseCore const & ellipse
) const {
    double r = ellipse.getDeterminantRadius() * getScale();
    Eigen::MatrixXd m = afw::math::shapelets::HermiteEvaluator::computeInnerProductMatrix(
        getOrder(), getOrder(), 1.0 / r, 1.0 / r
    );
    return m;
}

void ShapeletModelBasis::_integrate(lsst::ndarray::Array<Pixel, 1, 1> const & vector) const {
    afwShapelets::HermiteEvaluator shapeletEvaluator(_order);
    vector.deep() = 0.0;
    shapeletEvaluator.fillIntegration(vector);
    vector.deep() *= _scale * _scale;
}

void ShapeletModelBasis::_evaluate(
    lsst::ndarray::Array<Pixel, 2, 1> const & matrix,
    CONST_PTR(Footprint) const & footprint,
    lsst::afw::geom::Ellipse const & ellipse
) const {
    afw::geom::Ellipse fullEllipse(ellipse);
    fullEllipse.scale(_scale);
    afwShapelets::HermiteEvaluator shapeletEvaluator(_order);
    afw::geom::AffineTransform transform = fullEllipse.getGridTransform();
    Footprint::SpanList::const_iterator const spanEnd = footprint->getSpans().end();
    Footprint::SpanList::const_iterator spanIter = footprint->getSpans().begin();
    lsst::ndarray::Array<Pixel, 2, 1>::Iterator pixelIter = matrix.begin();
    for (; spanIter != spanEnd; ++spanIter) {
        afw::detection::Span const & span = **spanIter;
        for (int x = span.getX0(); x <= span.getX1(); ++x, ++pixelIter) {
            shapeletEvaluator.fillEvaluation(*pixelIter, transform(afw::geom::Point2D(x, span.getY())));
        }
    }
}

void ShapeletModelBasis::_evaluateRadialProfile(
    lsst::ndarray::Array<Pixel,2,1> const & profile,
    lsst::ndarray::Array<Pixel const,1,1> const & radii
) const {
    afwShapelets::BasisEvaluator basisEval(_order, afwShapelets::LAGUERRE);
    for (int r = 0; r < radii.getSize<0>(); ++r) {
        lsst::ndarray::Array<Pixel,1,1> row(profile[r]);
        basisEval.fillEvaluation(row, radii[r] / _scale, 0.0);
        for (int n = 0; n <= _order; ++n) {
            int offset = afwShapelets::computeOffset(n);
            for (int p = n, q = 0; q <= p; --p, ++q) {
                if (p != q) {
                    row[offset + 2 * q] = 0.0;
                    row[offset + 2 * q + 1] = 0.0;
                }
            }
        }
        afwShapelets::ConversionMatrix::convertOperationVector(
            row, afwShapelets::LAGUERRE, afwShapelets::HERMITE, _order
        );
    }
}

ModelBasis::Ptr ShapeletModelBasis::convolve(
    CONST_PTR(LocalPsf) const & psf
) const {
    if (psf->hasNativeShapelet()) {
        afwShapelets::MultiShapeletFunction s = psf->getNativeShapelet(afwShapelets::HERMITE);
        s.shiftInPlace(-afw::geom::Extent2D(psf->getPoint()));
        return convolve(s);
    } else {
        afwShapelets::ShapeletFunction s = 
            psf->computeShapelet(afwShapelets::HERMITE, getPsfShapeletOrder());
        s.getEllipse().getCenter() -= afw::geom::Extent2D(psf->getPoint());
        return convolve(s);
    }
}

ConvolvedShapeletModelBasis::Ptr ShapeletModelBasis::convolve(
    afwShapelets::ShapeletFunction const & psf
) const {
    return boost::make_shared<ConvolvedShapeletModelBasis>(*this, psf);
}

ModelBasis::Ptr ShapeletModelBasis::convolve(
    afwShapelets::MultiShapeletFunction const & psf
) const {
    CompoundShapeletModelBasis::ComponentVector components;
    components.push_back(boost::make_shared<ShapeletModelBasis>(_order, _scale));
    CompoundShapeletModelBasis::Ptr asCompound = CompoundShapeletBuilder(components).build();
    return asCompound->convolve(psf);
}

}}} // namespace lsst::meas::multifit
