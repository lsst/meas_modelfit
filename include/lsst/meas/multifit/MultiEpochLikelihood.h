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

#ifndef LSST_MEAS_MULTIFIT_MultiEpochLikelihood_h_INCLUDED
#define LSST_MEAS_MULTIFIT_MultiEpochLikelihood_h_INCLUDED

#include <vector>

#include "ndarray.h"

#include "lsst/pex/config.h"
#include "lsst/afw/geom/ellipses/Ellipse.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/shapelet/MultiShapeletBasis.h"

#include "lsst/meas/multifit/constants.h"
#include "lsst/meas/multifit/LogGaussian.h"
#include "lsst/meas/multifit/Likelihood.h"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief Control object used to initialize a MultiEpochObject.
 *
 *  Translated to Python as MultiEpochLikelihoodConfig; the Swig-wrapped C++ Control object can
 *  be created from the config object via the makeControl() method (see lsst.pex.config.wrap).
 */
class MultiEpochLikelihoodControl : public SingleEpochLikelihoodControl {
public:
    MultiEpochLikelihoodControl() : SingleEpochLikelihoodControl() {}
};

/**
 * An image at one epoch of a galaxy, plus associated info
 *
 * Includes one image of a galaxy and and associated footprint and multi-shapelet PSF model
 */
class EpochFootprint {
public:
    /**
     * @brief Construct a EpochFootprint
     *
     * @param[in] footprint     Footprint of source (galaxy) on calexp
     * @param[in] exposure      Subregion of calexp that includes footprint
     * @param[in] psfModel      Multi-shapelet representation of exposure PSF evaluated at location of galaxy
     */
    explicit EpochFootprint(
        lsst::afw::detection::Footprint const &footprint,
        lsst::afw::image::Exposure<Pixel> const &exposure,
        shapelet::MultiShapeletFunction const &psfModel
    );
    
    lsst::afw::detection::Footprint const footprint;  ///< footprint of source (galaxy)
    lsst::afw::image::Exposure<Pixel> const exposure; ///< subregion of exposure that includes footprint
    shapelet::MultiShapeletFunction const psfModel;   ///< multi-shapelet model of exposure PSF
    int const numPixels;                              ///< number of pixels in footprint
};

#ifndef SWIG
class EpochMatrixBuilder;   // defined and declared in MultiEpochLikelihood.cc
#endif

/// Likelihood class for use with multi-epoch modeling
class MultiEpochLikelihood : public Likelihood {
public:
    /// @copydoc Likelihood::getLinearDim
    virtual int getLinearDim() const { return _modelMatrix.getSize<1>(); }

    /// Return the sum of squares of the variance-weighted data vector
    virtual double getDataSquaredNorm() const { return _dataSquaredNorm; }

    /// @copydoc Likelihood::evaluate
    /// @note ellipse is in cooadd coordinates
    virtual LogGaussian evaluate(afw::geom::ellipses::Ellipse const & ellipse) const;

    /**
     * @brief Initialize the MultiEpochLikelihood
     *
     * @param[in] ctrl              Control object with various options
     * @param[in] basis             Basis object that defines the galaxy model to fit.
     * @param[in] coaddWcs          WCS of coadd
     * @param[in] sourceSkyPos      Sky position of source (galaxy)
     * @param[in] epochImageList    List of shared pointers to EpochFootprint
     */
    explicit MultiEpochLikelihood(
        MultiEpochLikelihoodControl const & ctrl,
        shapelet::MultiShapeletBasis const & basis,
        afw::image::Wcs const & coaddWcs,
        afw::coord::Coord const & sourceSkyPos,
        std::vector<PTR(EpochFootprint)> const & epochImageList
    );

private:
    int _totPixels;             ///< total number of pixels in all calexp
    double _dataSquaredNorm;    ///< sum of squares of _weightedData
    PixelArray1 _weights;       ///< vector of weights for all calexp concatenated
                                ///< = 1/sqrt(variance) pixels if ctrl.usePixelWeights
                                ///< = exp^mean(log(1/sqrt(variance))) if !ctrl.usePixelWeights,
                                ///<   where the mean is computed separately for each epoch
    PixelArray1 _weightedData;  ///< vector of weighted image pixels for all calexp concatenated
    PixelArray2CM _modelMatrix; ///< model matrix; set by evaluate
    std::vector<PTR(EpochMatrixBuilder)>  _epochMatrixBuilderList;
                                ///< list of shared_ptr to EpochMatrixBuilder;
                                ///< one entry for each element of epochImageList
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_MultiEpochLikelihood_h_INCLUDED
