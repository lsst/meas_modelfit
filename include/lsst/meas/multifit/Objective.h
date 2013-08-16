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

#ifndef LSST_MEAS_MULTIFIT_Objective_h_INCLUDED
#define LSST_MEAS_MULTIFIT_Objective_h_INCLUDED

#include "ndarray.h"

#include "lsst/pex/config.h"
#include "lsst/afw/geom/ellipses/Ellipse.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/shapelet/MultiShapeletBasis.h"

#include "lsst/meas/multifit/constants.h"
#include "lsst/meas/multifit/LogGaussian.h"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief Base class for optimizer/sampler objective functions that compute likelihood at a point.
 *
 *  Objective abstracts the problem of computing the likelihood over different kinds of data; we
 *  expect to have two major subclasses: one for single-epoch modeling, and one for multi-epoch
 *  modeling.
 */
class Objective
#ifndef SWIG
 : private boost::noncopyable
#endif
{
public:

    /// Return the number of linear dimensions
    virtual int getLinearDim() const = 0;

    /// Return the sum of squares of the variance-weighted data vector
    virtual double getDataSquaredNorm() const = 0;

    /**
     *  @brief Evaluate the likelihood at the given point given an ellipse.
     *
     *  Because we want to take advantage of the fact that the likelihood of the linear amplitude parameters
     *  is Gaussian, at each point in the nonlinear ellipse parameter space, we compute that Gaussian
     *  distribution in the amplitudes, in terms of the maximum likelihood vector, the Fisher matrix,
     *  and the sum of squared residuals at the maximum likelihood point.  See LogGaussian for more
     *  information.
     */
    virtual LogGaussian evaluate(afw::geom::ellipses::Ellipse const & ellipse) const = 0;

    virtual ~Objective() {}

};

/**
 *  @brief Control object used to initialize a SingleEpochObject.
 *
 *  Translated to Python as SingleEpochObjectiveConfig; the Swig-wrapped C++ Control object can
 *  be created from the config object via the makeControl() method (see lsst.pex.config.wrap).
 */
class SingleEpochObjectiveControl {
public:

    LSST_CONTROL_FIELD(usePixelWeights, bool,
                       "whether to individually weigh pixels using the variance image.");
    LSST_CONTROL_FIELD(useSVD, bool,
                       "Use SVD instead of Eigensystem decomposition to solve linear least-squares. "
                       "While both approaches are robust against rank-deficient matrices, and both "
                       "factorizations are performed in double precision, with the Eigensystem "
                       "method the normal equations will first be formed in single precision, which "
                       "makes them more subject to round-off error."
    );
    LSST_CONTROL_FIELD(useApproximateExp, bool,
                       "whether to use fast approximate exponentials when evaluating the model");

    SingleEpochObjectiveControl() : usePixelWeights(true), useSVD(false), useApproximateExp(false) {}
};

/// Objective class for use with single-epoch modeling
class SingleEpochObjective : public Objective {
public:

    /// @copydoc Objective::getLinearDim
    virtual int getLinearDim() const { return _modelMatrix.getSize<1>(); }

    /// Return the sum of squares of the variance-weighted data vector
    virtual double getDataSquaredNorm() const { return _dataSquaredNorm; }

    /// @copydoc Objective::evaluate
    virtual LogGaussian evaluate(afw::geom::ellipses::Ellipse const & ellipse) const;

    /**
     *  @brief Initialize the SingleEpochObjective
     *
     *  @param[in] ctrl      Control object with various options.
     *  @param[in] basis     Basis object that defines the galaxy model to fit.
     *  @param[in] psf       Multi-shapelet representation of the PSF evaluated at the
     *                       location of the galaxy.
     *  @param[in] image     MaskedImage to fit to.
     *  @param[in] footprint Footprint that defines the pixel region to include in the fit.
     */
    explicit SingleEpochObjective(
        SingleEpochObjectiveControl const & ctrl,
        shapelet::MultiShapeletBasis const & basis,
        shapelet::MultiShapeletFunction const & psf,
        afw::image::MaskedImage<Pixel> const & image,
        afw::detection::Footprint const & footprint
    );

private:
    double _dataSquaredNorm;
    PixelArray1 _weights;
    PixelArray1 _weightedData;
    PixelArray2CM _modelMatrix;
    shapelet::MultiShapeletMatrixBuilder<Pixel> _matrixBuilder;
};

namespace detail {
    #ifndef SWIG
    /**
     *  @brief A MultiShapeletMatrixBuilder factory with arguments suitable for SingleEpochObjective
     *
     *  @param[in] ctrl         SingleEpochObjective control object with various options.
     *  @param[in] basis        Basis object that defines the galaxy model to fit.
     *  @param[in] psf          Multi-shapelet representation of the PSF evaluated at the location of the galaxy.
     *  @param[in] footprint    Footprint that defines the pixel region to include in the fit.
     */
    shapelet::MultiShapeletMatrixBuilder<Pixel> makeShapeletMatrixBuilder(
        SingleEpochObjectiveControl const & ctrl,
        shapelet::MultiShapeletBasis const & basis,
        shapelet::MultiShapeletFunction const & psf,
        afw::detection::Footprint const & footprint
    );
    #endif
} // namespace ::detail

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_Objective_h_INCLUDED
