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

#ifndef LSST_MEAS_MULTIFIT_BaseEvaluator
#define LSST_MEAS_MULTIFIT_BaseEvaluator

#include "lsst/ndarray.h"
#include "lsst/meas/multifit/constants.h"
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief Represents a Bayesian prior on the coefficients given the parameters.
 *
 *  Priors on coefficients are represented as an unnormalized Gaussian factor multiplied
 *  by a completely arbitrary function of the coefficients:
 *  @f[
 *     P(\mu|\phi) = e^{-\frac{1}{2}(\mu-y)^T\!F(\mu-y)} g(\mu)
 *  @f]
 *  where @f$\mu@f$ are the coefficients, @f$y@f$ and @f$F@f$ are the Gaussian vector and
 *  matrix, and @f$g(\mu)@f$ is evaluated using operator()().
 *
 *  Note that the full prior must be normalized:
 *  @f[
 *     \int d\mu P(\mu|\phi) = 1
 *  @f]
 *  if the Bayesian evidence is to be computed correctly.  At the very least, the
 *  normalization must not vary with @f$\phi@f$.
 */
class CoefficientPrior : private boost::noncopyable {
public:

    typedef boost::shared_ptr<CoefficientPrior> Ptr;
    typedef boost::shared_ptr<CoefficientPrior const> ConstPtr;

    /// @brief Mean vector of the Gaussian factor of the prior (@f$y@f$) .
    Eigen::VectorXd const & getGaussianVector() const { return _gaussianVector; }

    /// @brief Inverse of the covariance matrix of the Gaussian factor of the prior (@f$F@f$).
    Eigen::MatrixXd const & getGaussianMatrix() const { return _gaussianMatrix; }

    /// @brief Evaluate the general factor of the prior (@f$g(\mu)@f$).
    virtual double operator()(lsst::ndarray::Array<Pixel const,1,1> const & coefficients) const = 0;

    virtual ~CoefficientPrior() {}

protected:
    Eigen::VectorXd _gaussianVector;
    Eigen::MatrixXd _gaussianMatrix;
};

/**
 *  @brief Code interface for nested linear models.
 *
 *  BaseEvaluator models an objective function of the form
 *  @f[
 *      q(\mu,\phi) = \frac{1}{2} (A\mu-x)^T\!(A\mu-x) 
 *           + \frac{1}{2}\ln\left|2\pi\Sigma^{-1}\right|
 *           + \frac{1}{2} (\mu-y)^T\!F(\mu-y)
 *           - \ln g(\mu)
 *  @f]
 *  where
 *   - @f$\mu@f$ are the linear parameters (called "coefficients") of the model
 *   - @f$\phi@f$ are the nonlinear parameters (called simply "parameters") of the model.
 *   - @f$A(\phi)@f$ is the model matrix that defines a linear model for the coefficients
 *     given the parameters.  If per-pixel weights are used, @f$A@f$ should include them
 *     (i.e. @f$A@f$ is the original model matrix with each row divided by the corresponding
 *     pixel's uncertainty).
 *   - @f$x@f$ is the data vector.  If per-pixel weights are used, @f$x@f$ includes them.
 *   - @f$\Sigma@f$ is the pixel uncertainty covariance matrix, assumed to be diagonal.
 *   - @f$y(\phi)@f$ is the mean of the Gaussian factor of a Bayesian prior on the coefficients.
 *     (see CoefficientPrior).
 *   - @f$F(\phi)@f$ is the inverse of the covariance matrix of the Gaussian factor of a 
 *     Bayesian prior on the coefficients (see CoefficientPrior).
 *   - @f$g(\mu,\phi)@f$ is the arbitrary factor in a Bayesian prior on the coefficents
 *     (see CoefficientPrior).
 *
 *  Note that @f$A@f$, @f$y@f$, @f$F@f$, and @f$g@f$ are implicitly functions of $\phi$
 *  in the above equation.
 *
 *  @f$q@f$ arises as the negative logarithm of the product of the likelihood and coefficient prior:
 *  @f[
 *      q(\mu,\phi) \equiv -\ln \left[ P(x|\mu,\phi) P(\mu|\phi) \right]
 *  @f]
 *  so
 *  @f[
 *      \int d\mu\,e^{-q(\mu,\phi)} = P(x|\phi)
 *  @f]
 *  The intent of the evaluator interface is thus to isolate the Gaussian or
 *  nearly Gaussian dependency of the likelihood on the coefficients, allowing this integral
 *  to be evaluated efficiently.
 *
 *  BaseEvaluator subclasses should be immutable.
 */
class BaseEvaluator {
public:

    typedef boost::shared_ptr<BaseEvaluator> Ptr;

    /// @brief Size of data vector (number of rows of model matrix).
    virtual int getPixelCount() const = 0;

    /// @brief Size of coefficient vector (number of colums of model matrix).
    virtual int getCoefficientCount() const = 0;

    /// @brief Number of parameters.
    virtual int getParameterCount() const = 0;

    /**
     *  @brief Return the natural log of the normalization term of the Gaussian likelihood.
     *
     *  This is a constant that does not depend on the parameters or coefficients, and only
     *  matters when the Bayesian evidence is computed.  For per-pixel uncertainties @f$\sigma_i@f$
     *  or, equivalently, a diagonal covariance matrix @f$\Sigma@f$, the returned value is
     *  @f[
     *    \sum_i \ln \frac{2\pi}{\sigma_i} = \frac{1}{2}\ln \left|2\pi\Sigma^{-1}\right|
     *  @f]
     */
    virtual double getLogPixelErrorSum() const = 0;

    /**
     *  @brief Data vector.
     *
     *  If the data vector is weighted (divided by sigma) the evaluted model matrix should be as well.
     */
    virtual lsst::ndarray::Array<Pixel const,1,1> getDataVector() const = 0;

    /// @brief Return true if all parameters are in-bounds.
    bool checkBounds(lsst::ndarray::Array<double const,1,1> & parameters) const;

    /**
     *  @brief Clip the given parameter vector to the valid range and return a penalty
     *         that scales with how far the parameter vector was beyond the constraints.
     */
    double clipToBounds(lsst::ndarray::Array<double,1,1> const & parameters) const;

    /**
     *  @brief Evaluate the parameter-dependent products of the evaluator.
     *
     *  @param[out] modelMatrix  The model matrix @f$A@f$; an array to fill with shape 
     *                           (getDataSize(), getCoefficientSize().
     *  @param[in]  parameters   The parameters @f$\phi@f$; an array of parameters with
     *                           size getParameterSize().
     *
     *  @returns The Bayesian prior on the coefficients given the parameters @f$P(\mu|\phi)@f$.  The
     *           returned object may be invalidated or updated by a subsequent call to evaluate().
     *
     *  If the data vector is weighted, the output matrix should be as well (each row should be divided
     *  by the corresponding pixel sigma value).
     */
    CoefficientPrior::ConstPtr evaluate(
        lsst::ndarray::Array<Pixel,2,2> const & modelMatrix,
        lsst::ndarray::Array<double const,1,1> const & parameters
    ) const;

    /// @brief Fill the given array with the initial parameter vector.
    void writeInitialParameters(lsst::ndarray::Array<double,1,1> const & parameters) const;

    /**
     *  @brief Compute the integral @f$\int d\mu\,e^{-q(\mu,\phi)} = P(x|\phi)@f$.
     *
     *  The integral is computed using importance sampling with a Gaussian distribution
     *  formed from the product of the likelihood and the Gaussian factor of the coefficient
     *  prior.
     *  
     *  @param[in]  engine        Random number generator.
     *  @param[out] coefficients  (samples)x(coefficient count) array to fill with Monte Carlo samples.
     *  @param[out] weights       Normalized weights (sum to one) corresponding to the Monte Carlo samples.
     *  @param[in]  parameters    Parameter vector to evaluate at.
     *
     *  @returns a Monte Carlo estimate of @f$P(x|\phi)@f$
     */
    double integrate(
        Random & engine,
        lsst::ndarray::Array<Pixel,2,2> const & coefficients,
        lsst::ndarray::Array<Pixel,1,1> const & weights,
        lsst::ndarray::Array<double const,1,1> const & parameters
    ) const;

    virtual ~BaseEvaluator() {}

protected:

    virtual CoefficientPrior::ConstPtr _evaluate(
        ndarray::Array<Pixel,2,2> const & matrix,
        ndarray::Array<double const,1,1> const & parameters
    ) const = 0;

    virtual void _writeInitialParameters(ndarray::Array<double,1,1> const & parameters) const = 0;

    virtual double _clipToBounds(lsst::ndarray::Array<double,1,1> const & parameters) const = 0;

    virtual bool _checkBounds(lsst::ndarray::Array<double const,1,1> & parameters) const = 0;

private:
    void operator=(BaseEvaluator const &) {}
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_BaseEvaluator
