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

#ifndef LSST_MEAS_MULTIFIT_priors_h_INCLUDED
#define LSST_MEAS_MULTIFIT_priors_h_INCLUDED

#include "lsst/meas/multifit/constants.h"
#include "lsst/meas/multifit/LogGaussian.h"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief Base class for Bayesian priors.
 *
 *  If @f$\alpha@f$ are the amplitudes, @f$\theta@f$ are the nonlinear parameters, and @f$D@f$ is
 *  the data, then this class represents @f$P(\alpha,\theta)@f$, by providing a function to evaluate
 *  @f[
 *    P(\theta|D)\,P(D) = \int\!P(D|\alpha,\theta)\,P(\alpha,\theta)\,d\alpha
 *  @f]
 *  at fixed @f$\theta@f$.  Because @f$\alpha@f$ are linear parameters, @f$P(D|\alpha,\theta)@f$
 *  is Gaussian in @f$\alpha@f$, and because @f$\theta@f$ is fixed, it's usually convenient to
 *  think of the integral as:
 *  @f[
 *    P(\theta|D)\,P(D) = P(\theta)\int\!P(D|\alpha,\theta)\,P(\alpha|\theta)\,d\alpha
 *      = P(\theta)\,P(D|\theta)\,P(D)
 *  @f]
 *  Thus, we marginalize the likelihood in @f$\alpha@f$ at fixed @f$\theta@f$, and then multiply
 *  by the prior on @f$\theta@f$.
 */
class Prior {
public:

    /**
     *  @brief Marginalize over the amplitude likelihood at the given point in nonlinear parameter space,
     *         and multiply the result by the prior on the nonlinear paramters.
     *
     *  @param[in]  likelihood     The likelihood @f$P(D|\alpha,\theta)@f$ as a (nonnormalized)
     *                             Gaussian in @f$\alpha@f$ at fixed @f$\theta@f$.
     *  @param[in]  parameters     The nonlinear parameters @f$\theta@f$, considered fixed for
     *                             this calculation.
     *
     *  @return The marginal likelihood multiplied by the prior on @f$\theta@f$:
     *  @f[
     *    P(D|\theta)P(\theta) =
     *       P(\theta)\int\!P(D|\alpha,\theta)\,P(\alpha|\theta)\,d\alpha
     *  @f]
     *
     *  Note that while the joint likelihood is passed in in logarithmic parameters, the return value
     *  should be not be the log likelihood.
     */
    virtual double apply(LogGaussian const & likelihood, Vector const & parameters) const = 0;

    /**
     *  @brief Compute the nonnormalized expectation value of the flux at a nonlinear parameter point.
     *
     *  This is the integral
     *  @f[
     *     \int\!\text{flux}(\alpha,\theta)\,P(D|\alpha,\theta)\,P(\alpha,\theta)\,d\alpha
     *  @f]
     *  For models in which each amplitude is the flux in a component, then
     *  @f$\text{flux}(\alpha,\theta)=|\alpha|_1@f$
     */
    virtual double computeFluxExpectation(
        LogGaussian const & likelihood, Vector const & parameters
    ) const = 0;

    /**
     *  @brief Compute the nonnormalized expectation value of the squared flux at a
     *         nonlinear parameter point.
     *
     *  This is the integral
     *  @f[
     *     \int\!\left[\text{flux}(\alpha,\theta)\right]^2\,P(D|\alpha,\theta)\,P(\alpha,\theta)\,d\alpha
     *  @f]
     *  For models in which each amplitude is the flux in a component, then
     *  @f$\text{flux}(\alpha,\theta)=|\alpha|_1@f$
     */
    virtual double computeSquaredFluxExpectation(
        LogGaussian const & likelihood, Vector const & parameters
    ) const = 0;

    /**
     *  @brief Compute the nonnormalized expectation value of the fraction of flux in each component
     *         at a nonlinear parameter point.
     *
     *  This is the integral
     *  @f[
     *     \int\!\text{fraction}(\alpha,\theta)\,P(D|\alpha,\theta)\,P(\alpha,\theta)\,d\alpha
     *  @f]
     *  For models in which each amplitude is the flux in a component, then
     *  @f$\text{fraction}(\alpha,\theta)=\frac{\alpha}{|\alpha|_1}@f$
     */
    virtual Vector computeFractionExpectation(
        LogGaussian const & likelihood, Vector const & parameters
    ) const = 0;

    virtual ~Prior() {}

};

/**
 *  @brief A dead-simple prior for 2-component models that has no dependence on nonlinear parameters, and
 *         declares that models must be either all one component or the other.
 *
 *  This is a nonnormalized prior; it cannot be used for computing the Bayesian evidence in a sense
 *  that would be useful for model selection tests.
 *
 *  For a two-element amplitude vector @f$\alpha@f$, the value of the prior is:
 *  @f[
 *    P(\alpha) = \begin{cases}
 *       \beta & \text{for $\alpha_0 > 0, \alpha_1 = 0$}\\
 *       1-\beta & \text{for $\alpha_0 = 0, \alpha_1 > 0$}\\
 *       0 & \text{otherwise}
 *    \end{cases}
 *  @f]
 */
class SingleComponentPrior : public Prior {
public:

    explicit SingleComponentPrior(double beta) : _beta(beta) {}

    /// Return the fraction of objects expected to be pure component 0
    double getBeta() const { return _beta; }

    virtual double apply(LogGaussian const & likelihood, Vector const & parameters) const;

    virtual double computeFluxExpectation(
        LogGaussian const & likelihood, Vector const & parameters
    ) const;

    virtual double computeSquaredFluxExpectation(
        LogGaussian const & likelihood, Vector const & parameters
    ) const;

    virtual Vector computeFractionExpectation(
        LogGaussian const & likelihood, Vector const & parameters
    ) const;

private:
    double _beta;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_priors_h_INCLUDED
