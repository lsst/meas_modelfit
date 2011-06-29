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

#include "lsst/meas/multifit/BaseCoefficientPrior.h"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief Code interface for nested linear models.
 *
 *  BaseEvaluator models an objective function of the form
 *  @f[
 *      q(\mu,\phi) = \frac{1}{2} (A\mu-x)^T\!(A\mu-x) 
 *           + \frac{1}{2}\ln\left|2\pi\Sigma^{-1}\right|
 *           - \ln P(\mu|\phi)
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
 *   - @f$P(\mu|\phi)@f$ is a Bayesian prior on the coefficents (see CoefficientPrior).
 *
 *  Note that @f$A@f$, is implicitly a function of $\phi$ in the above equation.
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
    bool checkBounds(lsst::ndarray::Array<double const,1,1> const & parameters) const;

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
    BaseCoefficientPrior::ConstPtr evaluate(
        lsst::ndarray::Array<Pixel,2,2> const & modelMatrix,
        lsst::ndarray::Array<double const,1,1> const & parameters
    ) const;

    /// @brief Fill the given array with the initial parameter vector.
    void writeInitialParameters(lsst::ndarray::Array<double,1,1> const & parameters) const;

    /**
     *  @brief Compute the Monte Carlo integral @f$\int d\mu\,P(x|\mu,\phi)\,P(\mu|\phi) = P(x|\phi)@f$.
     *
     *  This delegates to BaseCoefficientPrior::integrate and multiplies by the constant pixel 
     *  covariance factor.
     *  
     *  @param[in]  engine        Random number generator.
     *  @param[out] coefficients  (samples)x(coefficient count) array to fill with Monte Carlo samples.
     *  @param[out] weights       (samples) array of normalized weights (sum to one).
     *  @param[in]  parameters    (parameter count) array of parameters @f$\phi@f$.
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

    virtual BaseCoefficientPrior::ConstPtr _evaluate(
        ndarray::Array<Pixel,2,2> const & matrix,
        ndarray::Array<double const,1,1> const & parameters
    ) const = 0;

    virtual void _writeInitialParameters(ndarray::Array<double,1,1> const & parameters) const = 0;

    virtual double _clipToBounds(lsst::ndarray::Array<double,1,1> const & parameters) const = 0;

    virtual bool _checkBounds(lsst::ndarray::Array<double const,1,1> const & parameters) const = 0;

private:
    void operator=(BaseEvaluator const &) {}
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_BaseEvaluator
