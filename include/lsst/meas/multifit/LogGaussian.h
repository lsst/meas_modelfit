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
 *  L(\alpha) = \frac{1}{2}r + g^T\alpha + \frac{1}{2}\alpha^T F \alpha
 *  @f]
 *  where @f$r=2L(0.0)@f$,
 *  @f$g=\left.\frac{\partial L}{\partial\alpha}\right|_{0}@f$ is the gradient, and
 *  @f$F=\frac{\partial^2 L}{\partial\alpha^2}@f$ is the Fisher matrix.
 *
 *  Because the constant @f$r@f$ has the same value for all samples in a SampleSet,
 *  it is is not saved with each LogGaussian even though it is logically associated
 *  with it (see SampleSet::getDataSquaredNorm()).
 */
class LogGaussian {
public:

    samples::Vector grad;        ///< gradient of log likelihood at @f$\alpha=0@f$
    samples::Matrix fisher;      ///< Fisher matrix (inverse of the covariance matrix)

    /// Evaluate the negative log-likelihood at the given point.
    double operator()(double r, samples::Vector const & alpha) const {
        samples::Vector delta = alpha - grad;
        return 0.5*r + 0.5*delta.dot(fisher*delta);
    }

    /// Initialize from parameters
    LogGaussian(samples::Vector const & grad_, samples::Matrix const & fisher_) :
        grad(grad_), fisher(fisher_)
    {}

    /// Initialize with zeros of the given dimension
    explicit LogGaussian(int dim) :
        grad(samples::Vector::Zero(dim)), fisher(samples::Matrix::Zero(dim, dim))
    {}

};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_LogGaussian_h_INCLUDED
