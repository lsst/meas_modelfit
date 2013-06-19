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
 *  SampleSet::computeExpectation() computes an appoximation to the integral
 *  @f[
 *    \int\!\int\! f(\theta)\,P(\theta|D)\,d\theta
 *  @f]
 *  where @f$P(\theta|D)@f$ is the marginal posterior defined by the SampleSet, as
 *  @f[
 *    \frac{1}{P(D)N}\sum_{n=1}^N \frac{f(\cdot,\theta_n)m_n}{q_n}
 *  @f]
 *
 *  It is the responsibility of ExpectationFunctor simply to compute @f$f(\theta)@f$.
 */
class ExpectationFunctor {
public:

    /// Initialize the ExpectationFunctor with the dimensionality of its result.
    explicit ExpectationFunctor(int outputDim_) : outputDim(outputDim_) {}

    /// Compute the expectation quantity at the given sample point; see ExpectationFunctor.
    virtual samples::Vector operator()(samples::VectorCMap const & parameters) const = 0;

    virtual ~ExpectationFunctor() {}

    /// The dimensionality of the quantity we wish to compute the expectation of.
    int const outputDim;

};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_ExpectationFunctor_h_INCLUDED
