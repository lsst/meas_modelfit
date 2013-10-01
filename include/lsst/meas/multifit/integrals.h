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

#ifndef LSST_MEAS_MULTIFIT_integrals_h_INCLUDED
#define LSST_MEAS_MULTIFIT_integrals_h_INCLUDED

#include "Eigen/Core"
#include "lsst/meas/multifit/constants.h"

namespace lsst { namespace meas { namespace multifit {

namespace detail {

double phid(double z);
double bvnu(double h, double k, double rho);

} // namespace detail

/**
 *  @brief Compute a multidimensional Gaussian integral over the space of all nonnegative vectors,
 *         returning the log of the integral for stability.
 *
 *  This routine computes:
 *  @f[
 *  -\ln\int_0^{\inf} dx e^{-g^T x - \frac{1}{2}x^T F x}
 *  @f]
 *
 *  @param[in]  grad        The 'g' vector in the above formula; so-called because it's usually
 *                          the gradient vector of a log-likelihood function.
 *  @param[in]  fisher      The 'F' matrix in the above formula; so-called because it's usually
 *                          the Fisher (second-derivative) matrix of a log-likelihood function.
 *
 *  Currently only the 1-d and 2-d cases are supported.
 */
double integrateGaussian(Vector const & grad, Matrix const & fisher);

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_integrals_h_INCLUDED
