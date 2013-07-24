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

#ifndef LSST_MEAS_MULTIFIT_MultiEpochObjective_h_INCLUDED
#define LSST_MEAS_MULTIFIT_MultiEpochObjective_h_INCLUDED

#include <vector>

#include "ndarray.h"

#include "lsst/pex/config.h"
#include "lsst/afw/geom/ellipses/Ellipse.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/shapelet/MultiShapeletBasis.h"

#include "lsst/meas/multifit/constants.h"
#include "lsst/meas/multifit/LogGaussian.h"
#include "lsst/meas/multifit/Objective.h"

namespace lsst { namespace meas { namespace multifit {

/**
 * An image at one epoch of a galaxy, plus associated info
 *
 * Includes one image of a galaxy and and associated footprint and multi-shapelet PSF model
 */
class EpochImage {
public:
    /**
     * @brief Construct a EpochImage
     *
     * @param[in] footprint     Footprint of source (galaxy) on calexp
     * @param[in] exposure      Subregion of calexp that includes footprint
     * @param[in] psfModel      Multi-shapelet representation of exposure PSF evaluated at location of galaxy
     */
    explicit EpochImage(
        lsst::afw::detection::Footprint const &footprint,
        lsst::afw::image::Exposure<Pixel> const &exposure,
        shapelet::MultiShapeletFunction const &psfModel
    );
    
    lsst::afw::detection::Footprint footprint;  ///< footprint of source (galaxy)
    lsst::afw::image::Exposure<Pixel> exposure; ///< subregion of exposure that includes footprint
    shapelet::MultiShapeletFunction psfModel;   ///< multi-shapelet model of exposure PSF
    int numPixels;                              ///< number of pixels in footprint
};


#ifndef SWIG
/**
* Contains a MatrixBuilder for one EpochImage, plus additional information
*
* Intended to be constructed and used internally by MultiEpochObjective.
*/
class EpochMatrixBuilder {
public:
    /**
    * Construct a EpochMatrixBuilder
    *
    * @param[in] begIndex       Bginning index of this component in _weights, etc.
    * @param[in] numPixels      Number of pixels in footprint of this component
    * @param[in] coaddToCalexp  Affine transform of coadd pixels->calexp pixels
    * @param[in] matrixBuilder  Multi-shapelet matrix builder for this component
    * 
    */
    explicit EpochMatrixBuilder(
        int begIndex,
        int numPixels,
        lsst::afw::geom::AffineTransform coaddToCalexp,
        shapelet::MultiShapeletMatrixBuilder<Pixel> matrixBuilder
    );

    int begIndex;           ///< beginning index of this component in _weights, etc.
    int numPixels;          ///< number of pixels in the footprint for this component
    lsst::afw::geom::AffineTransform coaddToCalexp; /// affine transform of coadd pixels->calexp pixels
    shapelet::MultiShapeletMatrixBuilder<Pixel> matrixBuilder;  ///< multishapelet matrix builder
};
#endif

/// Objective class for use with multi-epoch modeling
class MultiEpochObjective : public Objective {
public:
    /// @copydoc Objective::getLinearDim
    virtual int getLinearDim() const { return _modelMatrix.getSize<1>(); }

    /// Return the sum of squares of the variance-weighted data vector
    virtual double getDataSquaredNorm() const { return _dataSquaredNorm; }

    /// @copydoc Objective::evaluate
    /// @note ellipse is in cooadd coordinates
    virtual LogGaussian evaluate(afw::geom::ellipses::Ellipse const & ellipse) const;

    /**
     * @brief Initialize the MultiEpochObjective
     *
     * @param[in] ctrl              Control object with various options
     * @param[in] basis             Basis object that defines the galaxy model to fit.
     * @param[in] coaddWcs          WCS of coadd
     * @param[in] sourceSkyPos      Sky position of source (galaxy)
     * @param[in] epochImageList    List of shared pointers to EpochImage
     */
    explicit MultiEpochObjective(
        SingleEpochObjectiveControl const & ctrl,
        shapelet::MultiShapeletBasis const & basis,
        afw::image::Wcs const & coaddWcs,
        afw::coord::Coord const & sourceSkyPos,
        std::vector<CONST_PTR(EpochImage)> const & epochImageList
    );

private:
    lsst::afw::coord::Coord _sourceSkyPos;  ///< sky position of source (likely not needed)
    int _totPixels;             ///< total number of pixels in all calexp
    double _dataSquaredNorm;    ///< sum of squares of _weightedData
    PixelArray1 _weights;       ///< vector of weights for all calexp concatenated
                                ///< = 1/sqrt(variance) pixels if ctrl.usePixelWeights
                                ///< = exp^mean(log(1/sqrt(variance))) if !ctrl.usePixelWeights,
                                ///<   where the mean is computed separately for each component
    PixelArray1 _weightedData;  ///< vector of weighted image pixels for all calexp concatenated
    PixelArray2CM _modelMatrix; ///< model matrix; set by evaluate
    std::vector<CONST_PTR(EpochMatrixBuilder)>  _epochMatrixBuilderList;
                                ///< list of shared_ptr to EpochMatrixBuilder;
                                ///< one entry for each element of epochImageList
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_MultiEpochObjective_h_INCLUDED
