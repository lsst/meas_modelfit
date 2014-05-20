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

#ifndef LSST_MEAS_MULTIFIT_UnitTransformedLikelihood_h_INCLUDED
#define LSST_MEAS_MULTIFIT_UnitTransformedLikelihood_h_INCLUDED

#include <vector>
#include "boost/scoped_ptr.hpp"

#include "ndarray.h"

#include "lsst/pex/config.h"
#include "lsst/afw/geom/ellipses/Ellipse.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/detection/Footprint.h"

#include "lsst/meas/multifit/common.h"
#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/Likelihood.h"
#include "lsst/meas/multifit/UnitSystem.h"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief Control object used to initialize a UnitTransformedLikelihood.
 *
 *  Translated to Python as UnitTransformedLikelihoodConfig; the Swig-wrapped C++ Control object can
 *  be created from the config object via the makeControl() method (see lsst.pex.config.wrap).
 */
class UnitTransformedLikelihoodControl {
public:

    LSST_CONTROL_FIELD(usePixelWeights, bool,
                       "whether to individually weigh pixels using the variance image.");
    LSST_CONTROL_FIELD(useApproximateExp, bool,
                       "whether to use fast approximate exponentials when evaluating the model");

    UnitTransformedLikelihoodControl() : usePixelWeights(true), useApproximateExp(false) {}

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
     * @param[in] psf           Multi-shapelet representation of exposure PSF evaluated at location of galaxy
     */
    explicit EpochFootprint(
        afw::detection::Footprint const &footprint,
        afw::image::Exposure<Pixel> const &exposure,
        shapelet::MultiShapeletFunction const &psf
    );

    afw::detection::Footprint const footprint;  ///< footprint of source (galaxy)
    afw::image::Exposure<Pixel> const exposure; ///< subregion of exposure that includes footprint
    shapelet::MultiShapeletFunction const psf;   ///< multi-shapelet model of exposure PSF
};

/// Likelihood class for use with multi-epoch modeling
class UnitTransformedLikelihood : public Likelihood {
public:

    /// @copydoc Likelihood::computeModelMatrix
    virtual void computeModelMatrix(
        ndarray::Array<Pixel,2,-1> const & modelMatrix,
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        bool doApplyWeights=true
    ) const;

    /**
     * @brief Initialize a UnitTransformedLikelihood with data from multiple exposures.
     *
     * @param[in] model             Object that defines the model to fit and its parameters.
     * @param[in] fixed             Model parameters that are held fixed.
     * @param[in] fitSys            Geometric and photometric system to fit in
     * @param[in] position          Sky position of object being fit
     * @param[in] epochFootprintList   List of shared pointers to EpochFootprint
     * @param[in] ctrl              Control object with various options
     */
    explicit UnitTransformedLikelihood(
        PTR(Model) model,
        ndarray::Array<Scalar const,1,1> const & fixed,
        UnitSystem const & fitSys,
        afw::coord::Coord const & position,
        std::vector<PTR(EpochFootprint)> const & epochFootprintList,
        UnitTransformedLikelihoodControl const & ctrl
    );

    /**
     * @brief Initialize a UnitTransformedLikelihood with data from multiple exposures.
     *
     * @param[in] model             Object that defines the model to fit and its parameters.
     * @param[in] fixed             Model parameters that are held fixed.
     * @param[in] fitSys            Geometric and photometric system to fit in
     * @param[in] position          Sky position of object being fit
     * @param[in] exposure          Exposure containing the data to fit
     * @param[in] footprint         Footprint that defines the pixels to include in the fit
     * @param[in] psf               Shapelet approximation to the PSF
     * @param[in] ctrl              Control object with various options
     */
    explicit UnitTransformedLikelihood(
        PTR(Model) model,
        ndarray::Array<Scalar const,1,1> const & fixed,
        UnitSystem const & fitSys,
        afw::coord::Coord const & position,
        afw::image::Exposure<Pixel> const & exposure,
        afw::detection::Footprint const & footprint,
        shapelet::MultiShapeletFunction const & psf,
        UnitTransformedLikelihoodControl const & ctrl
    );

    virtual ~UnitTransformedLikelihood();

private:
    class Impl;
    boost::scoped_ptr<Impl> _impl;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_UnitTransformedLikelihood_h_INCLUDED
