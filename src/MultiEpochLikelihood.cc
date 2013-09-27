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
#include "lsst/afw/image/Calib.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/detection/FootprintArray.cc"  // yes .cc; see the file for an explanation
#include "lsst/shapelet/MultiShapeletBasis.h"
#include "lsst/meas/multifit/MultiEpochLikelihood.h"

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
    afw::detection::Footprint const &footprint,
    afw::image::Exposure<Pixel> const &exposure,
    shapelet::MultiShapeletFunction const & psfModel
) :
    footprint(footprint),
    exposure(afw::image::Exposure<Pixel>(exposure, false)),
    psfModel(psfModel),
    numPixels(footprint.getNpix())
{}

MultiEpochLikelihood::MultiEpochLikelihood(
    PTR(Model) model,
    afw::image::Wcs const & coaddWcs,
    afw::image::Calib const & coaddCalib,
    afw::coord::Coord const & sourceSkyPos,
    std::vector<PTR(EpochFootprint)> const & epochImageList,
    MultiEpochLikelihoodControl const & ctrl
) : Likelihood(model) {
    int totPixels = std::accumulate(epochImageList.begin(), epochImageList.end(), 0, componentPixelSum);
    _data = ndarray::allocate(totPixels);
    afw::geom::AffineTransform coaddToSky = coaddWcs.linearizePixelToSky(sourceSkyPos, afw::geom::radians);
    int dataOffset = 0;
    for (
        std::vector<PTR(EpochFootprint)>::const_iterator imPtrIter = epochImageList.begin();
        imPtrIter != epochImageList.end(); ++imPtrIter
    ) {
        int const numPixels = (*imPtrIter)->numPixels; // used to shorten some expressions below

        afw::geom::AffineTransform skyToCalexp =
            (*imPtrIter)->exposure.getWcs()->linearizeSkyToPixel(sourceSkyPos, afw::geom::radians);

        _epochLikelihoods.push_back(
            boost::make_shared<SingleEpochLikelihood>(
                model, (**imPtrIter).psfModel,
                (**imPtrIter).exposure.getMaskedImage(),
                (**imPtrIter).footprint,
                skyToCalexp * coaddToSky,
                coaddCalib.getFluxMag0().first / (**imPtrIter).exposure.getCalib()->getFluxMag0().first,
                _data[ndarray::view(dataOffset, dataOffset + numPixels)],
                ctrl
            )
        );

        dataOffset += numPixels;
    }
}

void MultiEpochLikelihood::computeModelMatrix(
    ndarray::Array<Pixel,2,-1> const & modelMatrix,
    ndarray::Array<Scalar const,1,1> const & parameters
) const {
    int dataOffset = 0;
    for (
        EpochLikelihoodVector::const_iterator i = _epochLikelihoods.begin();
        i != _epochLikelihoods.end();
        ++i
    ) {
        (**i).computeModelMatrix(
            modelMatrix[ndarray::view(dataOffset, dataOffset + (**i).getDataDim())()],
            parameters
        );
        dataOffset += (**i).getDataDim();
    }
}

}}} // namespace lsst::meas::multifit
