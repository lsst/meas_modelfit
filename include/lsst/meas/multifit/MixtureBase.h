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

#ifndef LSST_MEAS_MULTIFIT_MixtureBase_h_INCLUDED
#define LSST_MEAS_MULTIFIT_MixtureBase_h_INCLUDED

#include "ndarray_fwd.h"

#include "lsst/base.h"
#include "lsst/afw/math/Random.h"
#include "lsst/afw/table/io/Persistable.h"
#include "lsst/meas/multifit/constants.h"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief Base class for Mixture probability distributions
 *
 *  This base class provides the common API used by AdaptiveImportanceSampler and SampleSet.
 *  Subclasses are specialized for the number of dimensions, and should be used directly when possible.
 */
class MixtureBase :
    public afw::table::io::PersistableFacade<MixtureBase>,
    public afw::table::io::Persistable {
public:

    /// Return the number of dimensions
    virtual int getDimension() const = 0;

    /// Return the number of dimensions
    virtual int getComponentCount() const = 0;

    /**
     *  @brief Evaluate the distribution probability density function (PDF) at the given points
     *
     *  @param[in] x       array of points, shape=(numSamples, N)
     *  @param[out] p      array of probability values, shape=(numSamples,)
     */
    virtual void evaluate(
        ndarray::Array<Scalar const,2,1> const & x,
        ndarray::Array<Scalar,1,0> const & p
    ) const = 0;

    /**
     *  @brief Evaluate the contributions of each component to the full probability at the given points
     *
     *  @param[in]  x     points to evaluate at, with number of columns equal to the number of dimensions
     *  @param[in]  p     array to fill, with number of columns equal to the number of components
     */
    virtual void evaluateComponents(
        ndarray::Array<Scalar const,2,1> const & x,
        ndarray::Array<Scalar,2,1> const & p
    ) const = 0;

    /**
     *  @brief Draw random variates from the distribution.
     *
     *  @param[in,out] rng random number generator
     *  @param[out] x      array of points, shape=(numSamples, N)
     */
    virtual void draw(afw::math::Random & rng, ndarray::Array<Scalar,2,1> const & x) const = 0;

    /**
     *  @brief Perform an Expectation-Maximization step, updating the component parameters to match
     *         the given weighted samples.
     *
     *  @param[in] x       array of variables, shape=(numSamples, N)
     *  @param[in] w       array of weights, shape=(numSamples,)
     *  @param[in] tau1    damping parameter (see below)
     *  @param[in] tau2    damping parameter (see below)
     *
     *  The updates to the @f$\sigma@f$ matrices are damped according to:
     *  @f[
     *  \sigma_d = \alpha\sigma_1 + (1-\alpha)\sigma_0
     *  @f]
     *  Where @f$\sigma_0@f$ is the previous matrix, @f$\sigma_1@f$ is the undamped update,
     *  and @f$\sigma_d@f$ is the damped update.  The parameter @f$\alpha@f$ is set
     *  by the ratio of the determinants:
     *  @f[
     *   r \equiv \frac{|\sigma_1|}{|\sigma_0|}
     *  @f]
     *  When @f$r \ge \tau_1@f$, @f$\alpha=1@f$; when @f$r \lt \tau_1@f$, it is rolled off
     *  quadratically to @f$\tau_2@f$.
     */
    virtual void updateEM(
        ndarray::Array<Scalar const,2,1> const & x,
        ndarray::Array<Scalar const,1,1> const & w,
        double tau1=0.0, double tau2=0.5
    ) = 0;

    /// Polymorphic deep copy
    virtual PTR(MixtureBase) clone() const = 0;

    virtual ~MixtureBase() {}

    inline friend std::ostream & operator<<(std::ostream & os, MixtureBase const & self) {
        self._stream(os);
        return os;
    }

protected:
    virtual void _stream(std::ostream & os) const = 0;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_MixtureBase_h_INCLUDED
