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

#ifndef LSST_MEAS_MULTIFIT_ExpectationFunctor_h_INCLUDED
#define LSST_MEAS_MULTIFIT_ExpectationFunctor_h_INCLUDED

#include "lsst/meas/multifit/SampleSet.h"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief Functor base class for use with computing expectation values of SampleSets.
 *
 *  SampleSet::computeExpectation() computes the integral
 *  @f[
 *    \int\!\int\! f(\alpha,\theta)\,P(\alpha,\theta|D)\,d\alpha\,d\theta
 *  @f]
 *  where @f$P(\alpha,\theta|D)@f$ is the joint posterior defined indirectly by the
 *  SampleSet and @f$f(\alpha,\theta)@f$ is the vector quantity whose expectation we want
 *  to compute, using the Monte Carlo approximation
 *  @f[
 *    \frac{1}{P(D)N}\sum_{n=1}^N \frac{1}{q_n}
 *            \int\!e^{-L_n(\alpha)}\,P(\alpha,\theta_n)\,f(\alpha,\theta_n)\,d\alpha
 *  @f]
 *  with
 *  @f[
 *    P(D) \approx \frac{1}{N}\sum_{n=1}^N \frac{m_n}{q_n}
 *  @f]
 *
 *  For expectation functions that depend on the amplitudes, this calculation is complex and
 *  involves the particular form of the prior @f$P(\alpha,\theta)@f$ on the amplitude; it is
 *  the job of the ExpectationFunctor to compute the integral
 *  @f[
 *    \int\!e^{-L_n(\alpha)}\,P(\alpha,\theta_n)\,f(\alpha,\theta_n)\,d\alpha
 *  @f]
 *  Because Prior is also a polymorphic class hierarchy, this requires some degree of
 *  double-dispatch between arbitrary Priors and arbitrary ExpectationFunctors.
 *
 *  However, for expectation functions that do not depend on the amplitudes, @f$f(\cdot,\theta)@f$
 *  can be brought outside the integral, which is then just @f$m_n@f$, and the computation reduces
 *  to
 *  @f[
 *    \frac{1}{P(D)N}\sum_{n=1}^N \frac{f(\cdot,\theta_n)m_n}{q_n}
 *  @f]
 */
class ExpectationFunctor {
public:

    /// Initialize the ExpectationFunctor with the dimensionality of its result.
    explicit ExpectationFunctor(int outputDim_) : outputDim(outputDim_) {}

    /// Compute the expectation quantity at the given sample point; see ExpectationFunctor.
    virtual Eigen::VectorXd operator()(SamplePoint const & sample, Prior const & prior) const = 0;

    virtual ~ExpectationFunctor() {}

    /// The dimensionality of the quantity we wish to compute the expectation of.
    int const outputDim;

};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_ExpectationFunctor_h_INCLUDED
