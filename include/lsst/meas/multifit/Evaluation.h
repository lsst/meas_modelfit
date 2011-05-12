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

#ifndef LSST_MEAS_MULTIFIT_Evaluation
#define LSST_MEAS_MULTIFIT_Evaluation

#include "lsst/meas/multifit/BaseEvaluator.h"
#include "lsst/meas/multifit/GaussianDistribution.h"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief A lazy-evaluation object to efficiently use a BaseEvaluator.
 *
 *  Evaluation is essentially a lazy solver for "separable" nonlinear least squares problems,
 *  in which the problem can be written as linear least squares problem with a parameterized
 *  matrix, extended to include an optional Bayesian prior.  More precisely, Evaluation
 *  supports problems of the form:
 *  @f[
 *      \min_{\phi,x} q = \frac{1}{2}(A(\phi)x - y)^T (A(\phi)x - y) 
 *          + \frac{1}{2}(x - \mu(\phi))^T \Sigma{-1}(\phi) (x - \mu(\phi))
 *          + \frac{1}{2}\ln\left|2\pi\Sigma(\phi)\right| - \ln P(\phi)
 *  @f]
 *  where:
 *   - @f$\phi@f$ is the vector of nonlinear parameters (hereafter simply called "parameters")
 *   - @f$x@f$ is the vector of linear parameters ("coefficients")
 *   - @f$y@f$ is the data vector
 *   - @f$A$@f is the model matrix
 *   - @f$\mu@f$ and @f$\Sigma@f$ are the mean and covariance of a Gaussian prior on
 *     the coefficients, conditional on the parameters: @f$P(x|\phi) ~ \mathcal{N}(\mu,\Sigma)@f$
 *   - @f$P(\phi)$ is the prior on the parameters
 *   - @f$q@f$ is the objective value to be minimized
 *
 *  Some intermediate and derivative products produced by the Evaluation include:
 *   - the residuals vector @f$r = Ax - y@f$
 *   - the coefficient Fisher matrix @f$F = A^T A + \Sigma^{-1}@f$
 *   - the lower-triangular Cholesky factor @f$L@f$ of the coefficient Fisher matrix @f$L L^T = F@f$
 *   - partial derivatives of the model matrix and residuals vector with respect to the parameters
 *
 *  If constructed without a prior, only the first term is evaluated.
 *
 *  If constructed with a prior, the main distribution must define the prior on the parameters
 *  and its nested distribution (if it has one) must be a GaussianDistribution on the coefficients.
 */
class Evaluation : private boost::noncopyable {
public:

    /// @brief Construct with no prior and use the evaluator's initial parameters.
    Evaluation(
        BaseEvaluator::Ptr const & evaluator, bool robustSolver=false
    );

    /// @brief Construct with a prior and use the evaluator's initial parameters.
    Evaluation(
        BaseEvaluator::Ptr const & evaluator,
        BaseDistribution const & prior,
        bool robustSolver=false
    );

    /// @brief Construct with no prior and the given parameter vector.
    Evaluation(
        BaseEvaluator::Ptr const & evaluator,
        lsst::ndarray::Array<double const,1,1> const & parameters,
        bool robustSolver=false
    );

    /// @brief Construct with a prior and the given parameter vector.
    Evaluation(
        BaseEvaluator::Ptr const & evaluator,
        lsst::ndarray::Array<double const,1,1> const & parameters,
        BaseDistribution const & prior,
        bool robustSolver=false
    );

    /// @brief Construct with no prior and the given parameter vector.
    Evaluation(
        BaseEvaluator::Ptr const & evaluator,
        Eigen::VectorXd const & parameters,
        bool robustSolver=false
    );

    /// @brief Construct with a prior and the given parameter vector.
    Evaluation(
        BaseEvaluator::Ptr const & evaluator,
        Eigen::VectorXd const & parameters,
        BaseDistribution const & prior,
        bool robustSolver=false
    );

    /// @brief Update the parameters @f$\phi@f$.
    void update(lsst::ndarray::Array<double const,1,1> const & parameters);

    /// @brief Update the parameters @f$\phi@f$.
    void update(Eigen::VectorXd const & parameters);

    /// @brief Update both the parameters @f$\phi@f$ and the coefficients @f$x@f$.
    void update(
        lsst::ndarray::Array<double const,1,1> const & parameters, 
        lsst::ndarray::Array<double const,1,1> const & coefficients
    );

    /// @brief Update both the parameters @f$\phi@f$ and the coefficients @f$x@f$.
    void update(
        Eigen::VectorXd const & parameters,
        Eigen::VectorXd const & coefficients
    );

    /// @brief Explicitly set the coefficients @f$x@f$.
    void setCoefficients(lsst::ndarray::Array<double const,1,1> const & coefficients);

    /// @brief Explicitly set the coefficients @f$x@f$.
    void setCoefficients(Eigen::VectorXd const & coefficients);
    
    /// @brief Solve for the coefficients @f$x@f$ that minimize @f$q@f$ given @f$\phi@f$.
    void solveCoefficients();

    /// @brief Return the evaluator that defines the model matrix and data vector.
    BaseEvaluator::Ptr getEvaluator() const { return _evaluator; }

    /// @brief Return the distribution that defines the prior.  May be empty.
    BaseDistribution::ConstPtr getPrior() const { return _prior; }
    
    /// @brief Return the parameters @f$\phi@f$.
    lsst::ndarray::Array<double const,1,1> getParameters() const { return _parameters; }

    /**
     *  @brief The model matrix @f$A@f$ that maps coefficients to model values.
     *
     *  The order of dimensions is {data, coefficients}.
     */
    lsst::ndarray::Array<double const,2,2> getModelMatrix() const {
        ensureModelMatrix();
        return _modelMatrix;
    }

    /**
     *  @brief The derivative of the model matrix with respect to the parameters, 
     *         @f$\frac{\partial A}{\partial \phi}@f$.
     *
     *  The order of dimensions is {parameters, data, coefficients}.
     */
    lsst::ndarray::Array<double const,3,3> getModelMatrixDerivative() const {
        ensureModelMatrixDerivative();
        return _modelMatrixDerivative;
    }

    /**
     *  @brief The current coefficient vector @f$x@f$.
     *
     *  If the coefficients have not been explicitly set or solved for since the
     *  last parameter change, they will be solved for.
     */
    lsst::ndarray::Array<double const,1,1> getCoefficients() const {
        ensureCoefficients();
        return _coefficients;
    }

    /// @brief The model vector @f$r = Ax@f$.
    lsst::ndarray::Array<double const,1,1> getModelVector() const {
        ensureModelVector();
        return _modelVector;
    }

    /// @brief The residuals vector @f$r = Ax - y@f$.
    lsst::ndarray::Array<double const,1,1> getResiduals() const {
        ensureResiduals();
        return _residuals;
    }

    /**
     *  @brief The partial derivative of the residuals vector with respect to the parameters,
     *          @f$\frac{\partial r}{\partial \phi}@f$.
     *
     *  The dimensions are ordered {data, parameters}.
     *
     *  Note that this is a partial derivative, not a full derivative (which would include
     *  an additional term @f$\frac{\partial r}{\partial x} \frac{\partial x}{\partial \phi}@f$
     *  if the coefficients are solved for).
     */
    lsst::ndarray::Array<double const,2,2> getResidualsJacobian() const {
        ensureResidualsJacobian();
        return _residualsJacobian;
    }

    /// @brief The coefficient Fisher matrix $F = A^T A + \Sigma^{-1}$.
    lsst::ndarray::Array<double const,2,2> getCoefficientFisherMatrix() const {
        ensureCoefficientFisherMatrix();
        return _coefficientFisherMatrix;
    }

    /// @brief The lower-triangular Cholesky factor @f$L@f$ of @f$F@f$, $L L^T = F = A^T A + \Sigma^{-1}$.
    lsst::ndarray::Array<double const,2,2> getCoefficientFisherFactor() const {
        ensureCoefficientFisherFactor();
        return _coefficientFisherFactor;
    }

    /// @brief Return the objective value @f$q@f$.
    double getObjectiveValue() const {
        ensureObjectiveValue();
        return _objectiveValue;
    }

    ~Evaluation();

private:

#ifndef SWIG
    class LinearSolver;
    class CholeskySolver;
    class EigenSolver;
#endif

    void ensureModelMatrix() const;
    void ensureModelMatrixDerivative() const;
    void ensureCoefficients() const;
    void ensureModelVector() const;
    void ensureResiduals() const;
    void ensureResidualsJacobian() const;
    void ensureCoefficientFisherMatrix() const;
    void ensureCoefficientFisherFactor() const;
    void ensureObjectiveValue() const;

    void initialize();
    void updateNestedPrior();

    mutable int _status;
    BaseEvaluator::Ptr _evaluator;
    BaseDistribution::ConstPtr _prior;
    GaussianDistribution::Ptr _nestedPrior;
    boost::scoped_ptr<LinearSolver> _solver;
    mutable double _objectiveValue;
    ndarray::Array<double,1,1> _parameters;
    mutable ndarray::Array<double,2,2> _modelMatrix;
    mutable ndarray::Array<double,3,3> _modelMatrixDerivative;
    mutable ndarray::Array<double,1,1> _coefficients;
    mutable ndarray::Array<double,2,2> _coefficientFisherMatrix;
    mutable ndarray::Array<double,2,2> _coefficientFisherFactor;
    mutable ndarray::Array<double,1,1> _residuals;
    mutable ndarray::Array<double,2,2> _residualsJacobian;
    mutable ndarray::Array<double,1,1> _modelVector;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_Evaluation
