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

#ifndef LSST_MEAS_MODELFIT_optimizer_h_INCLUDED
#define LSST_MEAS_MODELFIT_optimizer_h_INCLUDED

#include "ndarray.h"

#include "lsst/base.h"
#include "lsst/pex/config.h"
#include "lsst/afw/table/Schema.h"
#include "lsst/meas/modelfit/common.h"
#include "lsst/meas/modelfit/Model.h"

namespace lsst { namespace meas { namespace modelfit {

class Likelihood;
class Prior;
class Optimizer;

/**
 *  @brief Base class for objective functions for Optimizer
 */
class OptimizerObjective {
public:

    int const dataSize;
    int const parameterSize;

    /**
     *  Return a concrete Objective object built from a Likelihood and Prior.
     *
     *  Most fitting problems that can be formulated in terms of
     *  (multi-shapelet) Models, Likelihoods, and Priors can just use this
     *  Objective.  The returned Objective relies on numerical derivatives,
     *  however, so simple problems where analytic derivatives are easy to
     *  implement may merit a custom OptimizerObjective.
     */
    static std::shared_ptr<OptimizerObjective> makeFromLikelihood(
        std::shared_ptr<Likelihood> likelihood,
        std::shared_ptr<Prior> prior = std::shared_ptr<Prior>()
    );

    /**
     *  Base class constructor; must be called by all subclasses.
     */
    OptimizerObjective(int dataSize_, int parameterSize_) :
        dataSize(dataSize_), parameterSize(parameterSize_)
    {}

    /**
     *  Evaluate the Objective on a 1-d grid.
     *
     *  This delegates to computeResiduals, and hence need not be reimplemented
     *  by derived classes.  It is intended mostly for diagnostic purposes;
     *  when investigating the behavior of the optimizer, it is frequently
     *  valuable to plot its path against the gridded objective function.
     *
     *  @param[in]  parameters    An array with shape (gridSize, parameterSize),
     *                            specifying the parameter vector at points on
     *                            the grid.
     *  @param[out] output        Output array for objective values with shape
     *                            (gridSize).  Must be allocated, but need not
     *                            be initialized.
     *
     *  Frequently, the arguments to this function will be flattened views into
     *  higher dimensional arrays, allowing it to be used to construct N-d
     *  grids despite taking only a 1-d sequence of grid points.
     */
    void fillObjectiveValueGrid(
        ndarray::Array<Scalar const,2,1> const & parameters,
        ndarray::Array<Scalar,1,1> const & output
    ) const;

    /**
     *  Evaluate the residuals of the model for a given parameter vector.
     *
     *  @param[in]  parameters    An array of parameters with shape (parameterSize).
     *  @param[out] residuals     Output array that will contain (model - data) on
     *                            return.  Must be allocated to shape (dataSize),
     *                            but need not be initialized.
     */
    virtual void computeResiduals(
        ndarray::Array<Scalar const,1,1> const & parameters,
        ndarray::Array<Scalar,1,1> const & residuals
    ) const = 0;

    /**
     *  Evaluate analytic derivatives of the model or signal that they are not available.
     *
     *  Subclasses that provide analytic derivatives should implement this method and
     *  return true.  Subclasses that do not should return false to indicate to the optimizer
     *  that numeric derivatives must be used instead.  Subclasses that can only sometimes
     *  compute analytic derivatives are not supported.
     *
     *  @param[in]  parameters    An array of parameters with shape (parameterSize).
     *  @param[out] derivatives   Output array that will contain d(model - data)/d(parameters) on
     *                            return.  Must be allocated to shape (dataSize, parameterSize),
     *                            but need not be initialized.
     */
    virtual bool differentiateResiduals(
        ndarray::Array<Scalar const,1,1> const & parameters,
        ndarray::Array<Scalar,2,-2> const & derivatives
    ) const {
        return false;
    }


    /**
     *  Return true if the Objective has a Bayesian prior as well as a likelihood.
     *
     *  The default implementation returns false.
     */
    virtual bool hasPrior() const { return false; }

    /**
     *  Compute the value of the Bayesian prior for the given parameter vector.
     *
     *  The default implementation simply returns 1.0 (appropriate for an unnormalized constant prior,
     *  which is mathematically equivalent to no prior).
     */
    virtual Scalar computePrior(ndarray::Array<Scalar const,1,1> const & parameters) const { return 1.0; }

    /**
     *  Compute the first and second derivatives of the Bayesian prior with respect to the parameters.
     *
     *  The default implementation simply sets the output arrays to 0.0 (appropriate for an constant prior,
     *  which is mathematically equivalent to no prior).
     *
     *  @param[in]  parameters    An array of parameters with shape (parameterSize).
     *  @param[out] gradient      First derivative of the prior with respect to the parameters.
     *                            Must be allocated to shape (parameterSize), but need not be initialized.
     *  @param[out] hessian       Second derivative of the prior with respect to the parameters.
     *                            Must be allocated to shape (parameterSize, parameterSize), but need
     *                            not be initialized.  This is a symmetric matrix, and only the lower
     *                            triangle is used.
     */
    virtual void differentiatePrior(
        ndarray::Array<Scalar const,1,1> const & parameters,
        ndarray::Array<Scalar,1,1> const & gradient,
        ndarray::Array<Scalar,2,1> const & hessian
    ) const {
        gradient.deep() = 0.0;
        hessian.deep() = 0.0;
    }

    virtual ~OptimizerObjective() {}
};

/**
 *  @brief Configuration object for Optimizer
 *
 *  Many of these configuration options pertain to how the trust region is
 *  updated.  It's easiest to understand these together rather than separately.
 *  At each iteration, a quadratic model of the objective function is formed.  We can use this
 *  model to predict how we expect the objective function to behave over a step, and compare it
 *  to how the actual objective function behaves.  To do this, we'll use the ratio of the
 *  actual reduction in the objective function to the predicted reduction in the objective function,
 *  and call this @f$\rho@f$.  Then,
 *   - the step is accepted, and the parameters updated, when @f$\rho >@f$ @c stepAcceptThreshold.
 *   - if @f$\rho > @f$ @c trustRegionGrowReductionRatio and the length of the step is greater than
 *     @c trustRegionGrowStepFraction times the current trust region radius, the trust region radius
 *     will be multiplied by @c trustRegionGrowFactor.
 *   - if @c trustRegionShrinkMinReductionRatio @f$< \rho < @f$ @c trustRegionShrinkMaxReductionRatio,
 *     the trust region radius will be multiplied by @c trustRegionShrinkFactor.
 */
class OptimizerControl {
public:
    LSST_CONTROL_FIELD(
        noSR1Term, bool,
        "If true, ignore the SR1 update term in the Hessian, resulting in a Levenberg-Marquardt-like method"
    );

    LSST_CONTROL_FIELD(
        skipSR1UpdateThreshold, double,
        "Skip the SR1 update if |v||s| / (|v||s|) is less than this threshold"
    );

    LSST_CONTROL_FIELD(
        minTrustRadiusThreshold, double,
        "If the trust radius falls below this threshold, consider the algorithm converged"
    );

    LSST_CONTROL_FIELD(
        gradientThreshold, double,
        "If the maximum of the gradient falls below this threshold, consider the algorithm converged"
    );

    LSST_CONTROL_FIELD(
        numDiffRelStep, double,
        "relative step size used for numerical derivatives (added to other steps)"
    );

    LSST_CONTROL_FIELD(
        numDiffAbsStep, double,
        "absolute step size used for numerical derivatives (added to other steps)"
    );

    LSST_CONTROL_FIELD(
        numDiffTrustRadiusStep, double,
        "step size (in units of trust radius) used for numerical derivatives (added to relative step)"
    );

    LSST_CONTROL_FIELD(
        stepAcceptThreshold, double,
        "steps with reduction ratio greater than this are accepted"
    );

    LSST_CONTROL_FIELD(
        trustRegionInitialSize, double,
        "the initial trust region will be set to this value"
    );

    LSST_CONTROL_FIELD(
        trustRegionGrowReductionRatio, double,
        "steps with reduction radio greater than this may increase the trust radius"
    );

    LSST_CONTROL_FIELD(
        trustRegionGrowStepFraction, double,
        "steps with length this fraction of the trust radius may increase the trust radius"
    );

    LSST_CONTROL_FIELD(
        trustRegionGrowFactor, double,
        "when increase the trust region size, multiply the radius by this factor"
    );

    LSST_CONTROL_FIELD(
        trustRegionShrinkReductionRatio, double,
        "steps with reduction radio less than this will decrease the trust radius"
    );

    LSST_CONTROL_FIELD(
        trustRegionShrinkFactor, double,
        "when reducing the trust region size, multiply the radius by this factor"
    );

    LSST_CONTROL_FIELD(
        trustRegionSolverTolerance, double,
        "value passed as the tolerance to solveTrustRegion"
    );

    LSST_CONTROL_FIELD(
        maxInnerIterations, int,
        "maximum number of iterations (i.e. function evaluations and trust region subproblems) per step"
    );

    LSST_CONTROL_FIELD(
        maxOuterIterations, int,
        "maximum number of steps"
    );

    LSST_CONTROL_FIELD(
        doSaveIterations, bool,
        "whether to save all iterations for debugging purposes"
    );

    OptimizerControl() :
        noSR1Term(false), skipSR1UpdateThreshold(1E-8),
        minTrustRadiusThreshold(1E-5),
        gradientThreshold(1E-5),
        numDiffRelStep(0.0), numDiffAbsStep(0.0), numDiffTrustRadiusStep(0.1),
        stepAcceptThreshold(0.0),
        trustRegionInitialSize(1.0),
        trustRegionGrowReductionRatio(0.75),
        trustRegionGrowStepFraction(0.8),
        trustRegionGrowFactor(2.0),
        trustRegionShrinkReductionRatio(0.25),
        trustRegionShrinkFactor(1.0/3.0),
        trustRegionSolverTolerance(1E-8),
        maxInnerIterations(20),
        maxOuterIterations(500),
        doSaveIterations(false)
    {}
};

class OptimizerHistoryRecorder {
public:

    OptimizerHistoryRecorder(
        afw::table::Schema & schema,
        std::shared_ptr<Model> model,
        bool doRecordDerivatives
    );

    explicit OptimizerHistoryRecorder(afw::table::Schema const & schema);

    void apply(
        int outerIterCount,
        int innerIterCount,
        afw::table::BaseCatalog & history,
        Optimizer const & optimizer
    ) const;

    void unpackDerivatives(
        ndarray::Array<Scalar const,1,1> const & nested,
        Vector & gradient,
        Matrix & hessian
    ) const;

    void unpackDerivatives(
        afw::table::BaseRecord const & record,
        Vector & gradient,
        Matrix & hessian
    ) const;

    void unpackDerivatives(
        ndarray::Array<Scalar const,1,1> const & nested,
        ndarray::Array<Scalar,1,1> const & gradient,
        ndarray::Array<Scalar,2,2> const & hessian
    ) const;

    void unpackDerivatives(
        afw::table::BaseRecord const & record,
        ndarray::Array<Scalar,1,1> const & gradient,
        ndarray::Array<Scalar,2,2> const & hessian
    ) const;

    void fillObjectiveModelGrid(
        afw::table::BaseRecord const & record,
        ndarray::Array<Scalar const,2,1> const & parameters,
        ndarray::Array<Scalar,1,1> const & output
    ) const;

    afw::table::Key<int> outer;
    afw::table::Key<int> inner;
    afw::table::Key<int> state;
    ScalarKey objective;
    ScalarKey prior;
    ScalarKey trust;
    ArrayKey parameters;
    ArrayKey derivatives;
};

/**
 *  @brief A numerical optimizer customized for least-squares problems with Bayesian priors
 *
 *  The algorithm used by Optimizer combines the Gauss-Newton approach of approximating
 *  the second-derivative (Hessian) matrix as the inner product of the Jacobian of the residuals, while
 *  maintaining a matrix of corrections to this to account for large residuals, which is updated
 *  using a symmetric rank-1 (SR1) secant formula.  We assume the prior has analytic first and second
 *  derivatives, but use numerical derivatives to compute the Jacobian of the residuals at every
 *  step.  A trust region approach is used to ensure global convergence.
 *
 *  We consider the function @f$f(x)@f$ we wish to optimize to have two terms, which correspond to
 *  negative log likelihood (@f$\chi^2/2=\|r(x)|^2@f$, where @f$r(x)@f$ is the vector of residuals
 *  at @f$x@f$) and negative log prior @f$q(x)=-\ln P(x)@f$:
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
class Optimizer {
public:

    typedef OptimizerObjective Objective;
    typedef OptimizerControl Control;
    typedef OptimizerHistoryRecorder HistoryRecorder;

    enum StateFlags {
        CONVERGED_GRADZERO = 0x0001,
        CONVERGED_TR_SMALL = 0x0002,
        CONVERGED = CONVERGED_GRADZERO | CONVERGED_TR_SMALL,
        FAILED_MAX_INNER_ITERATIONS = 0x0010,
        FAILED_MAX_OUTER_ITERATIONS = 0x0020,
        FAILED_MAX_ITERATIONS = 0x0030,
        FAILED_EXCEPTION = 0x0040,
        FAILED_NAN = 0x0080,
        FAILED = FAILED_MAX_INNER_ITERATIONS | FAILED_MAX_OUTER_ITERATIONS | FAILED_EXCEPTION | FAILED_NAN,
        STATUS_STEP_REJECTED = 0x0100,
        STATUS_STEP_ACCEPTED = 0x0200,
        STATUS_STEP = STATUS_STEP_REJECTED | STATUS_STEP_ACCEPTED,
        STATUS_TR_UNCHANGED = 0x1000,
        STATUS_TR_DECREASED = 0x2000,
        STATUS_TR_INCREASED = 0x4000,
        STATUS_TR = STATUS_TR_UNCHANGED | STATUS_TR_DECREASED | STATUS_TR_INCREASED,
        STATUS = STATUS_STEP | STATUS_TR,
    };

    Optimizer(
        std::shared_ptr<Objective const> objective,
        ndarray::Array<Scalar const,1,1> const & parameters,
        Control const & ctrl
    );

    std::shared_ptr<Objective const> getObjective() const { return _objective; }

    Control const & getControl() const { return _ctrl; }

    bool step() { return _stepImpl(0); }

    bool step(HistoryRecorder const & recorder, afw::table::BaseCatalog & history) {
        return _stepImpl(0, &recorder, &history);
    }

    int run() { return _runImpl(); }

    int run(HistoryRecorder const & recorder, afw::table::BaseCatalog & history) {
        return _runImpl(&recorder, &history);
    }

    int getState() const { return _state; }

    Scalar getObjectiveValue() const { return _current.objectiveValue; }

    ndarray::Array<Scalar const,1,1> getParameters() const { return _current.parameters; }

    ndarray::Array<Scalar const,1,1> getResiduals() const { return _current.residuals; }

    ndarray::Array<Scalar const,1,1> getGradient() const { return _gradient; }

    ndarray::Array<Scalar const,2,2> getHessian() const { return _hessian; }

    /// Remove the symmetric-rank-1 secant term from the Hessian, making it just (J^T J)
    void removeSR1Term();

private:

    struct IterationData {
        Scalar objectiveValue;
        Scalar priorValue;
        ndarray::Array<Scalar,1,1> parameters;
        ndarray::Array<Scalar,1,1> residuals;

        IterationData(int dataSize, int parameterSize);

        void swap(IterationData & other);
    };

    friend class OptimizerHistoryRecorder;

    bool _stepImpl(
        int outerIterCount,
        HistoryRecorder const * recorder=NULL,
        afw::table::BaseCatalog * history=NULL
    );

    int _runImpl(HistoryRecorder const * recorder=NULL, afw::table::BaseCatalog * history=NULL);

    void _computeDerivatives();

    int _state;
    std::shared_ptr<Objective const> _objective;
    Control _ctrl;
    double _trustRadius;
    IterationData _current;
    IterationData _next;
    ndarray::Array<Scalar,1,1> _step;
    ndarray::Array<Scalar,1,1> _gradient;
    ndarray::Array<Scalar,2,2> _hessian;
    ndarray::Array<Scalar,2,-2> _residualDerivative;
    Matrix _sr1b;
    Vector _sr1v;
    Vector _sr1jtr;
};

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
void solveTrustRegion(
    ndarray::Array<Scalar,1,1> const & x,
    ndarray::Array<Scalar const,2,1> const & F, ndarray::Array<Scalar const,1,1> const & g,
    double r, double tolerance
);

}}} // namespace lsst::meas::modelfit

#endif // !LSST_MEAS_MODELFIT_optimizer_h_INCLUDED
