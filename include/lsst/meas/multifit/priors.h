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

#include "lsst/base.h"
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
 *
 *  @note All return values of prior integral member functions are negative logarithms, to avoid
 *        floating-point underflow when expoentiating large numbers.
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
     *  @return The marginal negative log likelihood multiplied by the prior on @f$\theta@f$
     *          (i.e. a nonnormalized log posterior):
     *  @f[
     *    P(D|\theta)P(\theta) =
     *       P(\theta)\int\!P(D|\alpha,\theta)\,P(\alpha|\theta)\,d\alpha
     *  @f]
     *
     *  Note that both the joint likelihood passed in and the marginal likelihood returned
     *  are the negative natural log of the probabilities densities.
     */
    virtual double apply(LogGaussian const & likelihood, Vector const & parameters) const = 0;

    /**
     *  @brief Compute the nonnormalized negative log expectation value of the flux at a
     *         nonlinear parameter point.
     *
     *  This is the integral
     *  @f[
     *     -\ln\int\!\text{flux}(\alpha,\theta)\,P(D|\alpha,\theta)\,P(\alpha,\theta)\,d\alpha
     *  @f]
     *  For models in which each amplitude is the flux in a component, then
     *  @f$\text{flux}(\alpha,\theta)=|\alpha|_1@f$
     */
    virtual double computeFluxExpectation(
        LogGaussian const & likelihood, Vector const & parameters
    ) const = 0;

    /**
     *  @brief Compute the nonnormalized negative log expectation value of the squared flux at a
     *         nonlinear parameter point.
     *
     *  This is the integral
     *  @f[
     *     -\ln\int\!\left[\text{flux}(\alpha,\theta)\right]^2\,P(D|\alpha,\theta)\,P(\alpha,\theta)\,d\alpha
     *  @f]
     *  For models in which each amplitude is the flux in a component, then
     *  @f$\text{flux}(\alpha,\theta)=|\alpha|_1@f$
     */
    virtual double computeSquaredFluxExpectation(
        LogGaussian const & likelihood, Vector const & parameters
    ) const = 0;

    /**
     *  @brief Compute the nonnormalized negative log expectation value of the fraction of flux
     *         in each component at a nonlinear parameter point.
     *
     *  This is the integral
     *  @f[
     *     -\ln\int\!\text{fraction}(\alpha,\theta)\,P(D|\alpha,\theta)\,P(\alpha,\theta)\,d\alpha
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
 *  @brief A nonnormalized flat prior that throws exceptions when the amplitude likelihood
 *         unconstrained in one or more directions.
 *
 *  @note FlatPrior is a singleton, accessed only via the get() static member function.
 */
class FlatPrior : public Prior, private boost::noncopyable {
public:

    static PTR(FlatPrior) get();

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
    FlatPrior() {}
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_priors_h_INCLUDED
