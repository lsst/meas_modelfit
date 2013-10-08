// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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

#include <limits>

#include "ndarray/eigen.h"

#include "lsst/afw/detection/FootprintArray.cc"  // yes .cc; see the file for an explanation
#include "lsst/meas/multifit/SingleEpochLikelihood.h"

namespace lsst { namespace meas { namespace multifit {

SingleEpochLikelihood::SingleEpochLikelihood(
    PTR(Model) model,
    ndarray::Array<Scalar const,1,1> const & fixed,
    shapelet::MultiShapeletFunction const & psf,
    afw::image::MaskedImage<Pixel> const & image,
    afw::detection::Footprint const & footprint,
    SingleEpochLikelihoodControl const & ctrl
) :
    Likelihood(model, fixed), _needsTransform(false),
    _weights(afw::detection::flattenArray(footprint, image.getVariance()->getArray(), image.getXY0())),
    _ellipses(model->makeEllipseVector())
{
    _init(psf, image, footprint, ctrl);
}

SingleEpochLikelihood::SingleEpochLikelihood(
    PTR(Model) model,
    ndarray::Array<Scalar const,1,1> const & fixed,
    shapelet::MultiShapeletFunction const & psf,
    afw::image::MaskedImage<Pixel> const & image,
    afw::detection::Footprint const & footprint,
    afw::geom::AffineTransform const & transform,
    double fluxScaling,
    ndarray::Array<Pixel,1,1> const & dataBuffer,
    SingleEpochLikelihoodControl const & ctrl
) :
    Likelihood(model, fixed), _needsTransform(true),
    _weights(afw::detection::flattenArray(footprint, image.getVariance()->getArray(), image.getXY0())),
    _ellipses(model->makeEllipseVector()),
    _transform(transform)
{
    _data = dataBuffer;
    _init(psf, image, footprint, ctrl);
    _weights.deep() *= fluxScaling;
}

void SingleEpochLikelihood::_init(
    shapelet::MultiShapeletFunction const & psf,
    afw::image::MaskedImage<Pixel> const & image,
    afw::detection::Footprint const & footprint,
    SingleEpochLikelihoodControl const & ctrl
) {
    LSST_ASSERT_EQUAL(
        _fixed.getSize<0>(), _model->getFixedDim(),
        "Fixed parameter vector size (%d) does not match Model fixed parameter dimensionality (%d)",
        pex::exceptions::LengthErrorException
    );
    {
        ndarray::Array<Pixel,1,1> x = ndarray::allocate(footprint.getArea());
        ndarray::Array<Pixel,1,1> y = ndarray::allocate(footprint.getArea());
        int n = 0;
        for (
            afw::detection::Footprint::SpanList::const_iterator i = footprint.getSpans().begin();
            i != footprint.getSpans().end();
            ++i
        ) {
            for (afw::geom::Span::Iterator j = (**i).begin(); j != (**i).end(); ++j, ++n) {
                x[n] = j->getX();
                y[n] = j->getY();
            }
        }
        for (
            Model::BasisVector::const_iterator k = getModel()->getBasisVector().begin();
            k != getModel()->getBasisVector().end();
            ++k
        ) {
            _matrixBuilders.push_back(
                shapelet::MultiShapeletMatrixBuilder<Pixel>(**k, psf, x, y, ctrl.useApproximateExp)
            );
        }
    }
    _data = afw::detection::flattenArray(footprint, image.getImage()->getArray(), image.getXY0());
    // Convert from variance to weights (1/sigma); this is actually the usual inverse-variance
    // weighting, because we implicitly square it later.
    _weights.asEigen<Eigen::ArrayXpr>() = _weights.asEigen<Eigen::ArrayXpr>().sqrt().inverse();
    if (!ctrl.usePixelWeights) {
        // We want a single number for the weights, so we use the geometric mean, as that
        // preserves the determinant of the (diagonal) pixel covariance matrix.
        _weights.asEigen<Eigen::ArrayXpr>().setConstant(
            std::exp(_weights.asEigen<Eigen::ArrayXpr>().log().mean())
        );
    }
    _data.asEigen<Eigen::ArrayXpr>() *= _weights.asEigen<Eigen::ArrayXpr>();
}

void SingleEpochLikelihood::computeModelMatrix(
    ndarray::Array<Pixel,2,-1> const & modelMatrix,
    ndarray::Array<Scalar const,1,1> const & parameters
) const {
    int c = 0;
    _model->writeEllipses(parameters.begin(), _fixed.begin(), _ellipses.begin());
    for (std::size_t i = 0; i < _ellipses.size(); ++i) {
        if (_needsTransform) {
            _ellipses[i].transform(_transform).inPlace();
        }
        int k = _matrixBuilders[i].getBasisSize();
        _matrixBuilders[i].build(modelMatrix[ndarray::view()(c, c+k)], _ellipses[i]);
        c += k;
    }
}

}}} // namespace lsst::meas::multifit
