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
#include "boost/make_shared.hpp"
#include "ndarray/eigen.h"

#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/detection/FootprintArray.cc"  // yes .cc; see the file for an explanation
#include "lsst/shapelet/MultiShapeletBasis.h"
#include "lsst/meas/multifit/MultiEpochObjective.h"

namespace lsst { namespace meas { namespace multifit {

namespace {

    /**
     *  @brief Function intended for use with std algorithms to compute the cumulative sum
     *  of the number of pixels in a sequence of components
     *
     * Example:
     * numPixels = std::accumulate(epochImageList.begin(), epochImageList.end(), 0, componentPixelSum);
     *
     * @param[in] partialNumPixels: accumulated number of pixels, so far
     * @param[in] epochImagePtr: shared pointer to EpochFootprint whose numPixels is to be added to the total
     * @return partialNumPixels + component.numPixels
     */
    int componentPixelSum(int partialNumPixels, CONST_PTR(EpochFootprint) const &epochImagePtr) {
        return partialNumPixels + epochImagePtr->numPixels;
    }
    
} // anonymous

EpochFootprint::EpochFootprint(
    lsst::afw::detection::Footprint const &footprint,
    lsst::afw::image::Exposure<Pixel> const &exposure,
    shapelet::MultiShapeletFunction const & psfModel
) :
    footprint(footprint),
    exposure(afw::image::Exposure<Pixel>(exposure, false)),
    psfModel(psfModel),
    numPixels(footprint.getNpix())
{ }

/**
* Contains a MatrixBuilder for one EpochFootprint, plus additional information
*
* Intended to be constructed and used internally by MultiEpochObjective.
*/
class EpochMatrixBuilder {
public:
    /**
    * Construct a EpochMatrixBuilder
    *
    * @param[in] begIndex       Beginning index of this epoch's data in _weights, etc.
    * @param[in] numPixels      Number of pixels in this epoch's footprint
    * @param[in] coaddToCalexp  Affine transform of coadd pixels->calexp pixels
    * @param[in] matrixBuilder  Multi-shapelet matrix builder for this epoch
    * 
    */
    explicit EpochMatrixBuilder(
        int begIndex,
        int numPixels,
        lsst::afw::geom::AffineTransform coaddToCalexp,
        shapelet::MultiShapeletMatrixBuilder<Pixel> matrixBuilder
    );

    int begIndex;           ///< beginning index of this epoch's data in _weights, etc.
    int numPixels;          ///< number of pixels in this epoch's footprint
    lsst::afw::geom::AffineTransform coaddToCalexp; /// affine transform of coadd pixels->calexp pixels
    shapelet::MultiShapeletMatrixBuilder<Pixel> matrixBuilder;  ///< multishapelet matrix builder
};

EpochMatrixBuilder::EpochMatrixBuilder(
    int begIndex,
    int numPixels,
    lsst::afw::geom::AffineTransform coaddToCalexp,
    shapelet::MultiShapeletMatrixBuilder<Pixel> matrixBuilder
) :
    begIndex(begIndex),
    numPixels(numPixels),
    coaddToCalexp(coaddToCalexp),
    matrixBuilder(matrixBuilder)
{ }

MultiEpochObjective::MultiEpochObjective(
    MultiEpochObjectiveControl const & ctrl,
    shapelet::MultiShapeletBasis const & basis,
    afw::image::Wcs const & coaddWcs,
    afw::coord::Coord const & sourceSkyPos,
    std::vector<CONST_PTR(EpochFootprint)> const & epochImageList
) :
    _totPixels(std::accumulate(epochImageList.begin(), epochImageList.end(), 0, componentPixelSum)),
    _dataSquaredNorm(0),
    _weights(ndarray::allocate(_totPixels)),
    _weightedData(ndarray::allocate(_totPixels)),
    _modelMatrix(ndarray::allocate(_totPixels, basis.getSize())),
    _epochMatrixBuilderList()
{
    afw::geom::AffineTransform coaddToSky = coaddWcs.linearizePixelToSky(sourceSkyPos, afw::geom::radians);
    int begIndex = 0;
    for (std::vector<CONST_PTR(EpochFootprint)>::const_iterator imPtrIter = epochImageList.begin();
        imPtrIter != epochImageList.end(); ++imPtrIter) {

        begIndex = begIndex;
        int const numPixels = (*imPtrIter)->numPixels; // used to shorten some expressions below
        
        afw::geom::AffineTransform skyToCalexp =
            (*imPtrIter)->exposure.getWcs()->linearizeSkyToPixel(sourceSkyPos, afw::geom::radians);
        _epochMatrixBuilderList.push_back(
            boost::make_shared<EpochMatrixBuilder>(
                begIndex,
                numPixels,
                skyToCalexp * coaddToSky,
                detail::makeShapeletMatrixBuilder(ctrl, basis,
                    (*imPtrIter)->psfModel, (*imPtrIter)->footprint)
        ));
        
        afw::image::MaskedImage<Pixel> maskedImage = (*imPtrIter)->exposure.getMaskedImage();
        _weights.asEigen().segment(begIndex, numPixels) =
            afw::detection::flattenArray(
                (*imPtrIter)->footprint,
                maskedImage.getVariance()->getArray(),
                maskedImage.getXY0()).asEigen();
            
        if (!ctrl.usePixelWeights) {
            // the weight for this component is the same for all pixels: e^mean(log(weight))
            _weights.asEigen().segment(begIndex, numPixels).setConstant(
                std::exp(_weights.asEigen<Eigen::ArrayXpr>().segment(begIndex, numPixels).log().mean()));
        }

        _weightedData.asEigen().segment(begIndex, numPixels) =
            afw::detection::flattenArray(
                (*imPtrIter)->footprint,
                maskedImage.getImage()->getArray(),
                maskedImage.getXY0()).asEigen() * _weights.asEigen().segment(begIndex, numPixels);

        begIndex += numPixels;
    }
    _dataSquaredNorm = _weightedData.asEigen().cast<double>().squaredNorm();
}

LogGaussian MultiEpochObjective::evaluate(afw::geom::ellipses::Ellipse const & ellipse) const {
    int const numCols = _modelMatrix.getSize<1>();
    for (std::vector<CONST_PTR(EpochMatrixBuilder)>::const_iterator mbPtrIter =
        _epochMatrixBuilderList.begin(); mbPtrIter != _epochMatrixBuilderList.end(); ++mbPtrIter) {
        afw::geom::ellipses::Ellipse localEllipse = ellipse.transform((*mbPtrIter)->coaddToCalexp);

        int const endIndex = (*mbPtrIter)->begIndex + (*mbPtrIter)->numPixels;
        (*mbPtrIter)->matrixBuilder.build(
            _modelMatrix[ndarray::view((*mbPtrIter)->begIndex, endIndex)()],
            localEllipse);
    }

    // copied directly from SingleEpochObjective::evaluate
    _modelMatrix.asEigen<Eigen::ArrayXpr>().colwise() *= _weights.asEigen<Eigen::ArrayXpr>();
    LogGaussian result(_modelMatrix.getSize<1>());
    // grad and fisher are the first and second derivatives of the log-likelihood for a zero
    // amplitude vector, and they're also the terms in the normal equations we'd solve for
    // the maximum likelihood solution
    result.grad = -_modelMatrix.asEigen().adjoint().cast<samples::Scalar>()
        * _weightedData.asEigen().cast<samples::Scalar>();
    result.fisher.selfadjointView<Eigen::Lower>().rankUpdate(
        _modelMatrix.asEigen().adjoint().cast<samples::Scalar>(), 1.0
    );
    result.fisher = result.fisher.selfadjointView<Eigen::Lower>();
    return result;
}

}}} // namespace lsst::meas::multifit
