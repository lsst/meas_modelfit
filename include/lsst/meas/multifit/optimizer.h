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

#ifndef LSST_MEAS_MULTIFIT_optimizer_h_INCLUDED
#define LSST_MEAS_MULTIFIT_optimizer_h_INCLUDED

#include "Eigen/Core"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief Solve a symmetric quadratic matrix equation with a ball constraint.
 *
 *  This computes a near-exact solution to the "trust region subproblem" necessary
 *  in trust-region-based nonlinear optimizers:
 *  @f[
 *   \min_x{\quad g^T x + \frac{1}{2}x^T F x}\quad\quad\quad \text{s.t.} ||x|| \le r
 *  @f]
 *
 *  The tolerance parameter sets how close to @f$r@f$ we require the norm of the
 *  solution to be when it lies on the constraint, as a fraction of @f$r@f$ itself.
 *
 *  This implementation is based on the algorithm described in Section 4.3 of
 *  "Nonlinear Optimization" by Nocedal and Wright.
 */
Eigen::VectorXd solveTrustRegion(
    Eigen::MatrixXd const & F, Eigen::VectorXd const & g, double r, double tolerance
);

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_optimizer_h_INCLUDED
