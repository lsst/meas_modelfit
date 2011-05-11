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
#include "lsst/afw/math/shapelets/detail/HermiteConvolution.h"
#include "lsst/ndarray/eigen.h"

namespace mf = lsst::meas::multifit;
namespace afwShapelets = lsst::afw::math::shapelets;

namespace lsst { namespace meas { namespace multifit {

namespace { 

static int const SHAPELET_DEFAULT_ORDER = 8; // TODO: make this a policy value

class ConvolvedShapeletModelBasis : public ModelBasis {
public:

    int getOrder() const { return _convolution->getColOrder(); }

    double getScale() const { return _scale; }

    explicit ConvolvedShapeletModelBasis(
        ShapeletModelBasis const & basis,
        afwShapelets::ShapeletFunction const & psf
    ) : ModelBasis(basis.getSize()),
        _convolution(boost::make_shared<afwShapelets::detail::HermiteConvolution>(basis.getOrder(), psf)),
        _frontBasis(ShapeletModelBasis::make(_convolution->getRowOrder(), 1.0)),
        _scale(basis.getScale())
    {}

protected:

    virtual void _integrate(lsst::ndarray::Array<Pixel, 1, 1> const & vector) const {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Cannot integrate convolved basis."
        );      
    }

    virtual void _evaluate(
        ndarray::Array<Pixel, 2, 1> const & matrix,
        CONST_PTR(Footprint) const & footprint,
        afw::geom::Ellipse const & ellipse
    ) const {
        ndarray::Array<double,2,2> frontMatrix(
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

private:

    afwShapelets::detail::HermiteConvolution::Ptr _convolution;
    ShapeletModelBasis::Ptr _frontBasis;
    double _scale;
};

} // anonymous

}}} // namespace lsst::meas::multifit

void mf::ShapeletModelBasis::_integrate(lsst::ndarray::Array<Pixel, 1, 1> const & vector) const {
    afwShapelets::detail::HermiteEvaluator shapeletEvaluator(_order);
    vector.deep() = 0.0;
    shapeletEvaluator.fillIntegration(vector);
    vector.deep() *= _scale * _scale;
}

void mf::ShapeletModelBasis::_evaluate(
    lsst::ndarray::Array<double, 2, 1> const & matrix,
    CONST_PTR(Footprint) const & footprint,
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
        for (int x = span.getX0(); x <= span.getX1(); ++x, ++pixelIter) {
            shapeletEvaluator.fillEvaluation(*pixelIter, transform(afw::geom::Point2D(x, span.getY())));
        }
    }
}

mf::ModelBasis::Ptr mf::ShapeletModelBasis::convolve(
    CONST_PTR(LocalPsf) const & psf
) const {
    if (psf->hasNativeShapelet()) {
        afwShapelets::MultiShapeletFunction s = psf->getNativeShapelet(afwShapelets::HERMITE);
        s.shiftInPlace(-afw::geom::Extent2D(psf->getPoint()));
        return convolve(s);
    } else {
        afwShapelets::ShapeletFunction s = 
            psf->computeShapelet(afwShapelets::HERMITE, SHAPELET_DEFAULT_ORDER);
        s.getEllipse().getCenter() -= afw::geom::Extent2D(psf->getPoint());
        return convolve(s);
    }
}

mf::ModelBasis::Ptr mf::ShapeletModelBasis::convolve(
    afwShapelets::ShapeletFunction const & psf
) const {
    return boost::make_shared<ConvolvedShapeletModelBasis>(*this, psf);
}

mf::ModelBasis::Ptr mf::ShapeletModelBasis::convolve(
    afwShapelets::MultiShapeletFunction const & psf
) const {
    CompoundShapeletModelBasis::ComponentVector components;
    components.push_back(boost::make_shared<ShapeletModelBasis>(_order, _scale));
    CompoundShapeletModelBasis::Ptr asCompound = CompoundShapeletBuilder(components).build();
    return asCompound->convolve(psf);
}
