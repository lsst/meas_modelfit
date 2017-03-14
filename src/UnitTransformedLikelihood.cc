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
#include <algorithm>
#include <limits>
#include <numeric>

#include "boost/format.hpp"
#include <memory>
#include "ndarray/eigen.h"

#include "lsst/afw/image/Calib.h"
#include "lsst/afw/detection/FootprintArray.cc"  // yes .cc; see the file for an explanation
#include "lsst/shapelet/MatrixBuilder.h"
#include "lsst/meas/modelfit/UnitTransformedLikelihood.h"

namespace lsst { namespace meas { namespace modelfit {

namespace {

typedef std::vector< shapelet::MatrixBuilder<Pixel> > BuilderVector;
typedef std::vector< shapelet::MatrixBuilderFactory<Pixel> > FactoryVector;

/*
 * Function intended for use with std algorithms to compute the cumulative sum
 * of the number of pixels in a sequence of EpochFootprints
 */
int componentPixelSum(int partialNumPixels, CONST_PTR(EpochFootprint) const &epochImagePtr) {
    return partialNumPixels + epochImagePtr->footprint.getArea();
}

/*
 * Return a vector of MatrixBuilders, with one for each MultiShapeletBasis in the input vector,
 * using the pixel region defined by the given Footprint and the given shapelet PSF approximation.
 *
 * basisVector - vector of MultiShapeletBasis objects; will produce one MatrixBuilder for each.
 * psf - MultiShapeletFunction representation of the PSF
 * footprint - Footprint that defines the region of pixels that will be used in the fit.
 */
BuilderVector makeMatrixBuilders(
    Model::BasisVector const & basisVector,
    shapelet::MultiShapeletFunction const & psf,
    afw::detection::Footprint const & footprint
) {
    BuilderVector builders;
    FactoryVector factories;
    builders.reserve(basisVector.size());
    factories.reserve(basisVector.size());
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
    int workspaceSize = 0;
    for (Model::BasisVector::const_iterator k = basisVector.begin(); k != basisVector.end(); ++k) {
        factories.push_back(shapelet::MatrixBuilderFactory<Pixel>(x, y, **k, psf));
        workspaceSize = std::max(workspaceSize, factories.back().computeWorkspace());
    }
    shapelet::MatrixBuilderWorkspace<Pixel> workspace(workspaceSize);
    for (FactoryVector::const_iterator i = factories.begin(); i != factories.end(); ++i) {
        shapelet::MatrixBuilderWorkspace<Pixel> wsCopy(workspace); // share workspace between builders
        builders.push_back((*i)(wsCopy));
    }
    return builders;
}

/*
 *  Flatten image and variance arrays from a MaskedImage using a footprint, and transform
 *  the variance into weights.
 *
 *  image - MaskedImage whose image and variance pixels should be used in the fit
 *  footprint - Footprint that defines the pixels to be included in the fit
 *  data - array to be filled with flattened values from the MaskedImage's image plane
 *  weights - array to be filled with flattened values computed from the MaskedImage's variance plane
 *  usePixelWeights - if true, weights will be per-pixel inverse sqrt(variance); if false, a constant
 *                    average value will be used
 */
void setupArrays(
    afw::image::MaskedImage<Pixel> const & image,
    afw::detection::Footprint const & footprint,
    ndarray::Array<Pixel,1,1> const & data,
    ndarray::Array<Pixel,1,1> const & variance,
    ndarray::Array<Pixel,1,1> const & weights,
    ndarray::Array<Pixel,1,1> const & unweightedData,
    bool usePixelWeights
) {
    afw::detection::flattenArray(footprint, image.getImage()->getArray(), data, image.getXY0());
    afw::detection::flattenArray(footprint, image.getVariance()->getArray(), variance, image.getXY0());
    unweightedData.deep() = data;
    // Convert from variance to weights (1/sigma); this is actually the usual inverse-variance
    // weighting, because we implicitly square it later.
    weights.asEigen<Eigen::ArrayXpr>() = variance.asEigen<Eigen::ArrayXpr>().sqrt().inverse();
    if (!usePixelWeights) {
        // If we're not using per-pixel weights, we need to use a constant non-unit weight instead,
        // which we compute as the geometric mean of the per-pixel weights.  The choice of geometric
        // mean preserves the determinant of the covariance matrix and makes it irrelevant whether
        // we average the variances or average the weights, but there's no real statistical
        // motivation for making the weights uniform (we do it to prevent model bias) and hence no
        // rigorous choice.
        weights.deep() = std::exp(weights.asEigen<Eigen::ArrayXpr>().log().sum() / weights.getSize<0>());
    }
    data.asEigen<Eigen::ArrayXpr>() *= weights.asEigen<Eigen::ArrayXpr>();
}

} // anonymous

EpochFootprint::EpochFootprint(
    afw::detection::Footprint const &footprint_,
    afw::image::Exposure<Pixel> const &exposure_,
    shapelet::MultiShapeletFunction const & psf_
) :
    footprint(footprint_),
    exposure(afw::image::Exposure<Pixel>(exposure_, false)),
    psf(psf_)
{}

class UnitTransformedLikelihood::Impl {
public:

    class Epoch {
    public:

        Epoch(int nPix_, LocalUnitTransform const & transform_, BuilderVector const & builders_) :
            nPix(nPix_), transform(transform_), builders(builders_) {}

        int nPix;
        LocalUnitTransform transform;
        BuilderVector builders;
    };

    Impl() : scratch(afw::geom::ellipses::Quadrupole(), afw::geom::Point2D()) {}

    std::vector<Epoch> epochs;
    Model::EllipseVector ellipses;
    afw::geom::ellipses::Ellipse scratch;
};

UnitTransformedLikelihood::UnitTransformedLikelihood(
    PTR(Model) model,
    ndarray::Array<Scalar const,1,1> const & fixed,
    UnitSystem const & fitSys,
    afw::coord::Coord const & position,
    std::vector<PTR(EpochFootprint)> const & epochFootprintList,
    UnitTransformedLikelihoodControl const & ctrl
) : Likelihood(model, fixed), _impl(new Impl()) {
    int totPixels = std::accumulate(epochFootprintList.begin(), epochFootprintList.end(),
                                    0, componentPixelSum);
    _data = ndarray::allocate(totPixels);
    _variance = ndarray::allocate(totPixels);
    _weights = ndarray::allocate(totPixels);
    _unweightedData = ndarray::allocate(totPixels);
    _impl->epochs.reserve(epochFootprintList.size());
    _impl->ellipses = model->makeEllipseVector();
    int dataOffset = 0;
    for (
        std::vector<PTR(EpochFootprint)>::const_iterator imPtrIter = epochFootprintList.begin();
        imPtrIter != epochFootprintList.end();
        ++imPtrIter
    ) {
        int nPix = (**imPtrIter).footprint.getArea();
        int dataEnd = dataOffset + nPix;
        _impl->epochs.push_back(
            Impl::Epoch(
                nPix, LocalUnitTransform(position, fitSys, (**imPtrIter).exposure),
                makeMatrixBuilders(model->getBasisVector(), (**imPtrIter).psf, (**imPtrIter).footprint)
            )
        );
        setupArrays(
            (**imPtrIter).exposure.getMaskedImage(),
            (**imPtrIter).footprint,
            _data[ndarray::view(dataOffset, dataEnd)],
            _variance[ndarray::view(dataOffset, dataEnd)],
            _weights[ndarray::view(dataOffset, dataEnd)],
            _unweightedData[ndarray::view(dataOffset, dataEnd)],
            ctrl.usePixelWeights
        );
    }
}

UnitTransformedLikelihood::UnitTransformedLikelihood(
    PTR(Model) model,
    ndarray::Array<Scalar const,1,1> const & fixed,
    UnitSystem const & fitSys,
    afw::coord::Coord const & position,
    afw::image::Exposure<Pixel> const & exposure,
    afw::detection::Footprint const & footprint,
    shapelet::MultiShapeletFunction const & psf,
    UnitTransformedLikelihoodControl const & ctrl
) : Likelihood(model, fixed), _impl(new Impl()) {
    int totPixels = footprint.getArea();
    _data = ndarray::allocate(totPixels);
    _variance = ndarray::allocate(totPixels);
    _weights = ndarray::allocate(totPixels);
    _unweightedData = ndarray::allocate(totPixels);
    _impl->ellipses = model->makeEllipseVector();
    _impl->epochs.push_back(
        Impl::Epoch(
            totPixels, LocalUnitTransform(position, fitSys, exposure),
            makeMatrixBuilders(model->getBasisVector(), psf, footprint)
        )
    );
    setupArrays(exposure.getMaskedImage(), footprint, _data, _variance, _weights, _unweightedData,
                ctrl.usePixelWeights);
}

UnitTransformedLikelihood::~UnitTransformedLikelihood() {}

void UnitTransformedLikelihood::computeModelMatrix(
    ndarray::Array<Pixel,2,-1> const & modelMatrix,
    ndarray::Array<Scalar const,1,1> const & nonlinear,
    bool doApplyWeights
) const {
    getModel()->writeEllipses(nonlinear.begin(), _fixed.begin(), _impl->ellipses.begin());
    int dataOffset = 0;
    modelMatrix.deep() = 0.0;
    for (
        std::vector<Impl::Epoch>::const_iterator i = _impl->epochs.begin();
        i != _impl->epochs.end();
        ++i
    ) {
        int dataEnd = dataOffset + i->nPix;
        int amplitudeOffset = 0;
        for (std::size_t j = 0; j < _impl->ellipses.size(); ++j) {
            _impl->scratch = _impl->ellipses[j].transform(i->transform.geometric);
            int amplitudeEnd = amplitudeOffset + i->builders[j].getBasisSize();
            i->builders[j](
                modelMatrix[ndarray::view(dataOffset, dataEnd)(amplitudeOffset, amplitudeEnd)],
                _impl->scratch
            );
            amplitudeOffset = amplitudeEnd;
        }
        modelMatrix[ndarray::view(dataOffset, dataEnd)()] *= i->transform.flux;
        dataOffset = dataEnd;
    }
    if (doApplyWeights) {
        modelMatrix.asEigen<Eigen::ArrayXpr>().colwise() *= _weights.asEigen<Eigen::ArrayXpr>();
    }
}

}}} // namespace lsst::meas::modelfit
