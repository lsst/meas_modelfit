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

#ifndef LSST_MEAS_MULTIFIT_SingleEpochLikelihood_h_INCLUDED
#define LSST_MEAS_MULTIFIT_SingleEpochLikelihood_h_INCLUDED

#include <vector>

#include "ndarray.h"

#include "lsst/pex/config.h"
#include "lsst/afw/geom/ellipses/Ellipse.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/shapelet/MultiShapeletBasis.h"

#include "lsst/meas/multifit/Likelihood.h"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief Control object used to initialize a SingleEpochLikelihood.
 *
 *  Translated to Python as SingleEpochLikelihoodConfig; the Swig-wrapped C++ Control object can
 *  be created from the config object via the makeControl() method (see lsst.pex.config.wrap).
 */
class SingleEpochLikelihoodControl {
public:

    LSST_CONTROL_FIELD(usePixelWeights, bool,
                       "whether to individually weigh pixels using the variance image.");
    LSST_CONTROL_FIELD(useApproximateExp, bool,
                       "whether to use fast approximate exponentials when evaluating the model");

    SingleEpochLikelihoodControl() : usePixelWeights(true), useApproximateExp(false) {}
};

/// Likelihood class for use with single-epoch modeling
class SingleEpochLikelihood : public Likelihood {
public:

    /// @copydoc Likelihood::computeModelMatrix
    virtual void computeModelMatrix(
        ndarray::Array<Pixel,2,-1> const & modelMatrix,
        ndarray::Array<Scalar const,1,1> const & parameters
    ) const;

    /**
     *  @brief Initialize the SingleEpochLikelihood
     *
     *  This overload assumes the data frame is the same as the parameter frame (i.e. no
     *  transform or flux scaling is necessary).
     *
     *  @param[in] model        Object that defines the model to fit and its parameters
     *  @param[in] psf          Multi-shapelet representation of the PSF evaluated at the
     *                          location of the source.
     *  @param[in] image        MaskedImage to fit to.
     *  @param[in] footprint    Footprint that defines the pixel region to include in the fit.
     *  @param[in] ctrl      Control object with various options.
     */
    explicit SingleEpochLikelihood(
        PTR(Model) model,
        shapelet::MultiShapeletFunction const & psf,
        afw::image::MaskedImage<Pixel> const & image,
        afw::detection::Footprint const & footprint,
        SingleEpochLikelihoodControl const & ctrl
    );

    /**
     *  @brief Constructor intended for use by MultiEpochObjective
     *
     *  @param[in] model        Object that defines the model to fit and its parameters
     *  @param[in] psf          Multi-shapelet representation of the PSF evaluated at the
     *                          location of the source.
     *  @param[in] image        MaskedImage to fit to.
     *  @param[in] footprint    Footprint that defines the pixel region to include in the fit.
     *  @param[in] transform    Transform that maps parameter frame to data frame.
     *  @param[in] fluxScaling  Factor to multiply model matrix by in order to match the data's
     *                          photometric scaling
     *  @param[in] dataBuffer   Memory to be used for the data array (but filled by the constructor).
     *  @param[in] ctrl      Control object with various options.
     */
    explicit SingleEpochLikelihood(
        PTR(Model) model,
        shapelet::MultiShapeletFunction const & psf,
        afw::image::MaskedImage<Pixel> const & image,
        afw::detection::Footprint const & footprint,
        afw::geom::AffineTransform const & transform,
        double fluxScaling,
        ndarray::Array<Pixel,1,1> const & dataBuffer,
        SingleEpochLikelihoodControl const & ctrl
    );

private:

    typedef std::vector< shapelet::MultiShapeletMatrixBuilder<Pixel> > MatrixBuilderVector;

    void _init(
        shapelet::MultiShapeletFunction const & psf,
        afw::image::MaskedImage<Pixel> const & image,
        afw::detection::Footprint const & footprint,
        SingleEpochLikelihoodControl const & ctrl
    );

    bool _needsTransform;
    PixelArray1 _weights;
    MatrixBuilderVector _matrixBuilders;
    mutable Model::EllipseVector _ellipses;
    afw::geom::AffineTransform _transform;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_SingleEpochLikelihood_h_INCLUDED
