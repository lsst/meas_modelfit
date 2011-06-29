// -*- LSST-C++ -*-
/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010, 2011 LSST Corporation.
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

#ifndef LSST_MEAS_MULTIFIT_BaseCoefficientPrior
#define LSST_MEAS_MULTIFIT_BaseCoefficientPrior

#include "lsst/ndarray.h"
#include "lsst/meas/multifit/constants.h"
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief Represents a Bayesian prior on the coefficients given the parameters.
 *
 *  A BaseCoefficientPrior may be constructed at a particular parameter vector, but it
 *  is not required to store that parameter vector or be able to return it to the
 *  user; this allows the same object to be reused when the prior is independent
 *  of the parameters.
 *
 *  Note that the prior must be normalized:
 *  @f[
 *     \int d\mu P(\mu|\phi) = 1
 *  @f]
 *  if the Bayesian evidence is to be computed correctly.  At the very least, the
 *  normalization must not vary with @f$\phi@f$.
 */
class BaseCoefficientPrior : private boost::noncopyable {
public:

    typedef boost::shared_ptr<BaseCoefficientPrior> Ptr;
    typedef boost::shared_ptr<BaseCoefficientPrior const> ConstPtr;

    /// @brief Evaluate the value of the prior for a given coefficient vector.
    virtual double operator()(lsst::ndarray::Array<Pixel const,1,1> const & coefficients) const = 0;

    /**
     *  @brief Compute the Monte Carlo integral @f$\int d\mu\,P(x|\mu,\phi)\,P(\mu|\phi) = P(x|\phi)@f$.
     *
     *  The full likelihood @f$P(x|\mu,\phi)@f$ has the form
     *  @f[
     *      \frac{1}{2}(A\mu - x)^T(A\mu - x)
     *  @f]
     *  An additional normalization factor involving the data covariance matrix must be manually included
     *  by the user.  Subclass implementations should consider the possibility that @f$A@f$ is
     *  not full rank.
     *  
     *  @param[in]  engine        Random number generator.
     *  @param[out] coefficients  (samples)x(coefficient count) array to fill with Monte Carlo samples.
     *  @param[out] weights       (samples) array of normalized weights (sum to one).
     *  @param[in]  modelMatrix   (pixel count)x(coefficient count) model matrix @f$A@f$.
     *  @param[in]  dataVector    (pixel count) array of pixel values @f$x@f$.
     *
     *  @returns a Monte Carlo estimate of @f$P(x|\phi)@f$
     */
    virtual double integrate(
        Random & engine,
        lsst::ndarray::Array<Pixel,2,2> const & coefficients,
        lsst::ndarray::Array<Pixel,1,1> const & weights,
        lsst::ndarray::Array<double const,2,2> const & modelMatrix,
        lsst::ndarray::Array<double const,1,1> const & dataVector
    ) const = 0;

    virtual ~BaseCoefficientPrior() {}

};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_BaseCoefficientPrior
