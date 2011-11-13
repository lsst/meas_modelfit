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
#include <boost/shared_ptr.hpp>

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief A lazy-evaluation object to efficiently use a BaseEvaluator.
 *
 *  Evaluation is essentially a lazy solver for "separable" nonlinear least squares problems,
 *  in which the problem can be written as linear least squares problem with a parameterized
 *  matrix.  More precisely, Evaluation supports problems of the form:
 *  @f[
 *      \min_{\phi,x} q = \frac{1}{2}(A(\phi)x - y)^T (A(\phi)x - y)
 *  @f]

 *  where:
 *   - @f$\phi@f$ is the vector of nonlinear parameters (hereafter simply called "parameters")
 *   - @f$x@f$ is the vector of linear parameters ("coefficients")
 *   - @f$y@f$ is the data vector
 *   - @f$A@f$ is the model matrix
 *   - @f$q@f$ is the objective value to be minimized
 *
 *  Some intermediate and derivative products produced by the Evaluation include:
 *   - the residuals vector @f$r = Ax - y@f$
 *   - the coefficient Fisher matrix @f$F = A^T A@f$
 *   - partial derivatives of the model matrix and residuals vector with respect to the parameters
 *
 *  If constructed without a prior, only the first term is evaluated.
 *
 *  If constructed with a prior, the main distribution must define the prior on the parameters
 *  and its nested distribution (if it has one) must be a GaussianDistribution on the coefficients.
 */
class Evaluation : private boost::noncopyable {
public:
    typedef boost::shared_ptr<Evaluation> Ptr;
    typedef boost::shared_ptr<Evaluation const> ConstPtr;

    /// @brief Construct with no prior and use the evaluator's initial parameters.
    Evaluation(BaseEvaluator::Ptr const & evaluator, double const svThreshold=1e-8);

    /// @brief Construct with no prior and the given parameter vector.
    Evaluation(
        BaseEvaluator::Ptr const & evaluator,
        lsst::ndarray::Array<double const,1,1> const & parameters,
        double const svThreshold=1e-8
    );

    /// @brief Construct with no prior and the given parameter vector.
    Evaluation(
        BaseEvaluator::Ptr const & evaluator,
        Eigen::VectorXd const & parameters,
        double const svThreshold=1e-8
    );

    /// @brief Update the parameters @f$\phi@f$.
    void update(lsst::ndarray::Array<double const,1,1> const & parameters);

    /// @brief Update the parameters @f$\phi@f$.
    void update(Eigen::VectorXd const & parameters);

    /// @brief Update both the parameters @f$\phi@f$ and the coefficients @f$x@f$.
    void update(
        lsst::ndarray::Array<double const,1,1> const & parameters, 
        lsst::ndarray::Array<Pixel const,1,1> const & coefficients
    );

    /// @brief Update both the parameters @f$\phi@f$ and the coefficients @f$x@f$.
    void update(
        Eigen::VectorXd const & parameters,
        Eigen::VectorXd const & coefficients
    );

    /// @brief Explicitly set the coefficients @f$x@f$.
    void setCoefficients(lsst::ndarray::Array<Pixel const,1,1> const & coefficients);

    /// @brief Explicitly set the coefficients @f$x@f$.
    void setCoefficients(Eigen::VectorXd const & coefficients);
    
    /// @brief Solve for the coefficients @f$x@f$ that minimize @f$q@f$ given @f$\phi@f$.
    void solveCoefficients();

    /// @brief Return the evaluator that defines the model matrix and data vector.
    BaseEvaluator::Ptr getEvaluator() const { return _evaluator; }
    
    /// @brief Return the parameters @f$\phi@f$.
    lsst::ndarray::Array<double const,1,1> getParameters() const { return _parameters; }

    /**
     *  @brief The model matrix @f$A@f$ that maps coefficients to model values.
     *
     *  The order of dimensions is {data, coefficients}.
     */
    lsst::ndarray::Array<Pixel const,2,2> getModelMatrix() const {
        ensureModelMatrix();
        return _modelMatrix;
    }

    /**
     *  @brief The derivative of the model matrix with respect to the parameters, 
     *         @f$\frac{\partial A}{\partial \phi}@f$.
     *
     *  The order of dimensions is {parameters, data, coefficients}.
     */
    lsst::ndarray::Array<Pixel const,3,3> getModelMatrixDerivative() const {
        ensureModelMatrixDerivative();
        return _modelMatrixDerivative;
    }

    /**
     *  @brief The current coefficient vector @f$x@f$.
     *
     *  If the coefficients have not been explicitly set or solved for since the
     *  last parameter change, they will be solved for.
     */
    lsst::ndarray::Array<Pixel const,1,1> getCoefficients() const {
        ensureCoefficients();
        return _coefficients;
    }

    /// @brief The model vector @f$r = Ax@f$.
    lsst::ndarray::Array<Pixel const,1,1> getModelVector() const {
        ensureModelVector();
        return _modelVector;
    }

    /// @brief The residuals vector @f$r = Ax - y@f$.
    lsst::ndarray::Array<Pixel const,1,1> getResiduals() const {
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
    lsst::ndarray::Array<Pixel const,2,2> getResidualsJacobian() const {
        ensureResidualsJacobian();
        return _residualsJacobian;
    }

    /// @brief The coefficient Fisher matrix @f$F = A^T A@f$.
    lsst::ndarray::Array<Pixel const,2,2> getCoefficientFisherMatrix() const {
        ensureCoefficientFisherMatrix();
        return _coefficientFisherMatrix;
    }

    /// @brief The coefficient Covariance matrix @f$\Sigma = F^+ = (A^T A)^+@f$.
    lsst::ndarray::Array<Pixel const,2,2> getCoefficientCovarianceMatrix() const {
        ensureCoefficientCovarianceMatrix();
        return _coefficientCovarianceMatrix;
    }

    /// @brief Return the objective value @f$q@f$.
    double getObjectiveValue() const {
        ensureObjectiveValue();
        return _objectiveValue;
    }
    
    ~Evaluation();

private:

    void ensureModelMatrix() const;
    void ensureModelMatrixDerivative() const;
    void ensureCoefficients() const;
    void ensureModelVector() const;
    void ensureResiduals() const;
    void ensureResidualsJacobian() const;
    void ensureCoefficientFisherMatrix() const;
    void ensureCoefficientCovarianceMatrix() const;
    void ensureObjectiveValue() const;
    void ensureFactorization() const;

    void initialize();

    class Factorization;

    mutable int _status;
    mutable int _products;
    BaseEvaluator::Ptr _evaluator;
    mutable boost::scoped_ptr<Factorization> _factorization;
    mutable double _objectiveValue;
    ndarray::Array<double,1,1> _parameters;
    mutable ndarray::Array<Pixel,2,2> _modelMatrix;
    mutable ndarray::Array<Pixel,3,3> _modelMatrixDerivative;
    mutable ndarray::Array<Pixel,1,1> _coefficients;
    mutable ndarray::Array<Pixel,2,2> _coefficientFisherMatrix;
    mutable ndarray::Array<Pixel,2,2> _coefficientCovarianceMatrix;
    mutable ndarray::Array<Pixel,1,1> _residuals;
    mutable ndarray::Array<Pixel,2,2> _residualsJacobian;
    mutable ndarray::Array<Pixel,1,1> _modelVector;
    double _svThreshold;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_Evaluation
