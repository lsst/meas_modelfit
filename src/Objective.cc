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

#include "lsst/afw/detection/FootprintArray.cc"
#include "lsst/shapelet/MultiShapeletBasis.h"
#include "lsst/meas/multifit/Objective.h"

namespace lsst { namespace meas { namespace multifit {

namespace {

// We need to build the ctor args for MultiShapeletMatrixBuilder from SingleEpochObjective ctor args,
// and we can't do that directly in the SingleEpochObjective ctor.
shapelet::MultiShapeletMatrixBuilder<Pixel> makeShapeletMatrixBuilder(
    SingleEpochObjectiveControl const & ctrl,
    shapelet::MultiShapeletBasis const & basis,
    shapelet::MultiShapeletFunction const & psf,
    afw::detection::Footprint const & footprint
) {
    PixelArray1 x = ndarray::allocate(footprint.getArea());
    PixelArray1 y = ndarray::allocate(footprint.getArea());
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
    return shapelet::MultiShapeletMatrixBuilder<Pixel>(basis, psf, x, y, ctrl.useApproximateExp);
}

} // anonymous

SingleEpochObjective::SingleEpochObjective(
    SingleEpochObjectiveControl const & ctrl,
    shapelet::MultiShapeletBasis const & basis,
    shapelet::MultiShapeletFunction const & psf,
    afw::image::MaskedImage<Pixel> const & image,
    afw::detection::Footprint const & footprint
) :
    _weights(afw::detection::flattenArray(footprint, image.getVariance()->getArray(), image.getXY0())),
    _weightedData(afw::detection::flattenArray(footprint, image.getImage()->getArray(), image.getXY0())),
    _modelMatrix(ndarray::allocate(footprint.getArea(), basis.getSize())),
    _leastSquares(
        ctrl.useSVD ? afw::math::LeastSquares::DIRECT_SVD : afw::math::LeastSquares::NORMAL_EIGENSYSTEM,
        basis.getSize()
    ),
    _matrixBuilder(makeShapeletMatrixBuilder(ctrl, basis, psf, footprint))
{
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
    _weightedData.asEigen<Eigen::ArrayXpr>() *= _weights.asEigen<Eigen::ArrayXpr>();
}

LogGaussian SingleEpochObjective::evaluate(afw::geom::ellipses::Ellipse const & ellipse) const {
    _matrixBuilder.build(_modelMatrix, ellipse);
    _modelMatrix.asEigen<Eigen::ArrayXpr>().colwise() *= _weights.asEigen<Eigen::ArrayXpr>();
    LogGaussian result(_leastSquares.getDimension());
    _leastSquares.setDesignMatrix(_modelMatrix, _weightedData);
    _leastSquares.setThreshold(std::numeric_limits<Pixel>::epsilon());
    result.mu = _leastSquares.getSolution().asEigen();
    result.fisher = _leastSquares.getFisherMatrix().asEigen();
    samples::Vector residuals = _modelMatrix.asEigen().cast<samples::Scalar>() * result.mu
        - _weightedData.asEigen().cast<samples::Scalar>();
    result.r = 0.5 * residuals.squaredNorm();
    return result;
}

}}} // namespace lsst::meas::multifit
