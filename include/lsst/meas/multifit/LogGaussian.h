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

#ifndef LSST_MEAS_MULTIFIT_LogGaussian_h_INCLUDED
#define LSST_MEAS_MULTIFIT_LogGaussian_h_INCLUDED

#include "Eigen/Core"

#include "lsst/meas/multifit/constants.h"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief A non-normalized Gaussian probability distribution represented as its negative log,
 *         for use in dealing with negative log likelihoods of linear parameters.
 *
 *  For a log likelihood @f$L(\alpha) = -\ln P(D|\alpha)@f$ in
 *  linear parameters @f$\alpha@f$, we solve for the maximum likelihood point @f$\mu@f$,
 *  so we can expand to second-order with no first-order term:
 *  @f[
 *  L(\alpha) = r
 *     + \frac{1}{2}\left(\alpha - \mu\right)^T F\left(\alpha - \mu\right)
 *  @f]
 *  where @f$r=L(\mu)@f$ and
 *  @f$F=\left.\frac{\partial^2 L}{\partial\alpha^2}\right|_{\mu}@f$ is
 *  the Fisher matrix.
 */
class LogGaussian {
public:

    Pixel r;           ///< negative log-likelihood at mu; @f$\chi^2/2@f$
    Vector mu;          ///< maximum likelihood point
    Matrix fisher;      ///< Fisher matrix at mu (inverse of the covariance matrix)

    /// Evaluate the negative log-likelihood at the given point.
    double operator()(Vector const & alpha) const {
        Vector delta = alpha - mu;
        return r + 0.5*delta.dot(fisher*delta);
    }

    /// Initialize from parameters
    LogGaussian(Pixel r_, Vector const & mu_, Matrix const & fisher_) : r(r_), mu(mu_), fisher(fisher_) {}


    /// Initialize with zeros of the given dimension
    LogGaussian(int dim) : r(0.0), mu(Vector::Zero(dim)), fisher(Matrix::Zero(dim, dim)) {}

};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_LogGaussian_h_INCLUDED
