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

#include "lsst/pex/config.h"

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

/**
 *  @brief Base class for objective functions for PosteriorOptimizer
 *
 */
class PosteriorOptimizerObjective {
public:

    int const parameterSize;
    int const dataSize;

    PosteriorOptimizerObjective(int parameterSize_, int dataSize_) :
        parameterSize(parameterSize_), dataSize(dataSize_)
    {}

    virtual void computeResiduals(Eigen::VectorXd const & parameters, Eigen::VectorXd & residuals) const = 0;

    virtual bool hasPrior() const { return false; }

    virtual double computePrior(Eigen::VectorXd const & parameters) const;

    virtual void differentiatePrior(
        Eigen::VectorXd const & parameters,
        Eigen::VectorXd & gradient,
        Eigen::MatrixXd & hessian
    ) const;

    virtual ~PosteriorOptimizerObjective() {}
};

/// Configuration object for PosteriorOptimizer
class PosteriorOptimizerControl {
public:
    LSST_CONTROL_FIELD(
        noSR1Term, bool,
        "If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method"
    );

    LSST_CONTROL_FIELD(
        numDiffRelStep, double,
        "relative step size used for numerical derivatives (added to absolute step)"
    );

    LSST_CONTROL_FIELD(
        numDiffAbsStep, double,
        "absolute step size used for numerical derivatives (added to relative step)"
    );

    PosteriorOptimizerControl() : noSR1Term(false), numDiffRelStep(1E-8), numDiffAbsStep(1E-8) {}
};

/**
 *  @brief A numerical optimizer customized for least-squares problems with Bayesian priors
 *
 *  The algorithm used by PosteriorOptimizer combines the Gauss-Newton approach of approximating
 *  the second-derivative (Hessian) matrix as the inner product of the Jacobian of the residuals, while
 *  maintaining a matrix of corrections to this to account for large residuals, which is updated
 *  using a symmetric rank-1 (SR1) secant formula.  We assume the prior has analytic first and second
 *  derivatives, but use numerical derivatives to compute the Jacobian of the residuals at every
 *  step.  A trust region approach is used to ensure global convergence.
 *
 *  We consider the function @f$f(x)@f$we wish to optimize to have two terms, which correspond to
 *  negative log likelihood (@f$\chi^2/2@f$) and negative log prior @f$q(x)=-\ln P(x)@f$:
 *  @f[
 *   f(x) = \frac{1}{2}\|r(x)\|^2 + q(x)
 *  @f]
 *  At each iteration @f$k@f$, we expand @f$f(x)@f$ in a Taylor series in @f$s=x_{k+1}-x_{k}@f$:
 *  @f[
 *   f(x) \approx m(s) = f(x_k) + g_k^T s + \frac{1}{2}s^T H_k s
 *  @f]
 *  where
 *  @f[
 *   g_k \equiv \left.\frac{\partial f}{\partial x}\right|_{x_k} = J_k^T r_k + \nabla q_k;\quad\quad
 *   J_k \equiv \left.\frac{\partial r}{\partial x}\right|_{x_k}
 *  @f]
 *  @f[
 *   H_k = J_k^T J_k + \nabla^2 q_k + B_k
 *  @f]
 *  Here, @f$B_k@f$ is the SR1 approximation term to the second derivative term:
 *  @f[
 *    B_k \approx \sum_i \frac{\partial^2 r^{(i)}_k}{\partial x^2}r^{(i)}_k
 *  @f]
 *  which we initialize to zero and then update with the following formula:
 *  @f[
 *   B_{k+1} = B_{k} + \frac{v v^T}{v^T s};\quad\quad v\equiv J^T_{k+1} r_{k+1} - J^T_k r_k
 *  @f]
 *  Unlike the more common rank-2 BFGS update formula, SR1 updates are not guaranteed to produce a
 *  positive definite Hessian.  This can result in more accurate approximations of the Hessian
 *  (and hence more accurate covariance matrices), but it rules out line-search methods and the simple
 *  dog-leg approach to the trust region problem.  As a result, we should require fewer steps to
 *  converge, but spend more time computing each step; this is ideal when we expect the time spent
 *  in function evaluation to dominate the time per step anyway.
 */
class PosteriorOptimizer {
public:

    typedef PosteriorOptimizerObjective Objective;
    typedef PosteriorOptimizerControl Control;

    PosteriorOptimizer(
        PTR(Objective const) objective,
        Eigen::VectorXd const & parameters,
        Control const & ctrl
    );

    PTR(Objective const) getObjective() const { return _objective; }

    Control const & getControl() const { return _control; }

    bool step();

    int getState() const;

    Eigen::VectorXd const & getResiduals() const;

    Eigen::VectorXd const & getGradient() const;

    Eigen::MatrixXd const & getHessian() const;

private:
    PTR(Objective const) _objective;
    Control _ctrl;
    Eigen::VectorXd _parameters;
    Eigen::VectorXd _testParameters;
    Eigen::VectorXd _residuals;
    Eigen::VectorXd _testResiduals;
    Eigen::MatrixXd _jacobian;
    Eigen::MatrixXd _testJacobian;
    Eigen::MatrixXd _gradient;
    Eigen::MatrixXd _hessian;
    Eigen::MatrixXd _sr1b;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_optimizer_h_INCLUDED
