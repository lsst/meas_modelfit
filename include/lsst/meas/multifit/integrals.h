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
#include "ndarray.h"
#include "lsst/afw/math/Random.h"
#include "lsst/meas/multifit/common.h"

namespace lsst { namespace meas { namespace multifit { namespace detail {

/**
 *  @brief Compute univariate normal probabilities
 *
 *  This function computes
 *  @f[
 *    \frac{1}{2\pi}\int_z^{\infty}dx\;e^{-x^2/(2z^2)}
 *  @f]
 *
 *  This is just a simple wrapper around boost::math::erfc, used to provide the particular form
 *  expected by bvnu.
 */
double phid(double z);

/**
 *  @brief Compute bivariate normal probabilities
 *
 *  This function computes
 *  @f[
 *    \frac{1}{2\pi\sqrt{1-\rho^2}}\int_h^{\infty}dx\int_k^{\infty}dy
 *        \;e^{-(x^2 - 2\rho x y + y^2)/(2(1-\rho^2))}
 *  @f]
 *
 *  It is a reimplementation of the "bvnu" MatLab routine by Alan Genz.  The original implementation
 *  can be found at http://www.math.wsu.edu/faculty/genz/homepage.  It is based on the algorithm described
 *  in:
 *
 *  *Drezner & Wesolowsky (1989), "On the computation of the bivariate normal integral",
 *     Journal of Statist. Comput. Simul. 35, pp. 101-107.*
 *
 *  A copy of Genz's FORTRAN routine of the same is included in the tests directory, as it has been used
 *  to generate the test reference data there.  It should generally not need to be compiled by users.
 *
 *  Most of the time, this function should be called via the methods in the TruncatedGaussian class;
 *  it is exposed publically only for testing purposes.
 */
double bvnu(double h, double k, double rho);

}}}} // namespace lsst::meas::multifit::detail

#endif // !LSST_MEAS_MULTIFIT_integrals_h_INCLUDED
