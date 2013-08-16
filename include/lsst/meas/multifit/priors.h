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
#include "lsst/meas/multifit/parameters.h"
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
class Prior : private boost::noncopyable {
public:

    /// Return the definition of the parameter vector the prior assumes
    ParameterDefinition const & getParameterDefinition() const { return *_parameterDefinition; }

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
     *  are the negative natural log of the probability densities.
     */
    virtual samples::Scalar apply(LogGaussian const & likelihood,
                                  samples::Vector const & parameters) const = 0;

    virtual ~Prior() {}

protected:
    explicit Prior(ParameterDefinition const & parameterDefinition) :
        _parameterDefinition(&parameterDefinition) {}
private:
    ParameterDefinition const * _parameterDefinition;
};

/**
 *  @brief A flat prior with upper bounds on radius and ellipticity (and implicit lower bounds),
 *         and a nonnegative requirement on component amplitudes.
 */
class FlatPrior : public Prior {
public:

    virtual samples::Scalar apply(LogGaussian const & likelihood, samples::Vector const & parameters) const;

    /**
     *  @brief Construct a FlatPrior.
     *
     *  @param[in] maxRadius   Maximum radius in pixel coordinates.
     *  @param[in] maxEllipticity   Maximum ellipticity (in ReducedShear parametrization).
     */
    explicit FlatPrior(double maxRadius, double maxEllipticity=1.0);

private:
    double _maxRadius;
    double _maxEllipticity;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_priors_h_INCLUDED
