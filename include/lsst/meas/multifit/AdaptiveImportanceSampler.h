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

#ifndef LSST_MEAS_MULTIFIT_AdaptiveImportanceSampler_h_INCLUDED
#define LSST_MEAS_MULTIFIT_AdaptiveImportanceSampler_h_INCLUDED

#include <map>

#include "lsst/meas/multifit/BaseSampler.h"
#include "lsst/meas/multifit/Mixture.h"
#include "lsst/meas/multifit/priors.h"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief Control object for one iteration of adaptive importance sampling
 */
class ImportanceSamplerControl {
public:
    LSST_CONTROL_FIELD(nSamples, int, "Number of Monte Carlo samples to draw");
    LSST_CONTROL_FIELD(nUpdateSteps, int, "Number of Expectation-Maximization update iterations");
    LSST_CONTROL_FIELD(tau1, double, "Damping parameter for E-M update (see MixtureBase::updateEM)");
    LSST_CONTROL_FIELD(tau2, double, "Damping parameter for E-M update (see MixtureBase::updateEM)");
    LSST_CONTROL_FIELD(
        targetPerplexity, double,
        "Minimum value for normalized perplexity after this iteration; if the actual value is less "
        "than this, this iteration will be repeated up to maxRepeat times until the target is met. "
        "In addition, if any previous iteration meets this target, this iteration will be skipped."
    );
    LSST_CONTROL_FIELD(
        maxRepeat, int,
        "Maximum number of times this iteration will be repeated to meet the perplexityTarget"
    );

    ImportanceSamplerControl() :
        nSamples(2000), nUpdateSteps(2), tau1(1E-4), tau2(0.5), targetPerplexity(1.0), maxRepeat(0)
    {}
};

/**
 *  @brief Sampler class that performs Monte Carlo sampling, while iteratively updating the
 *         analytic distribution from which points are drawn.
 *
 *  Between the iterations defined in the control object, the prior is applied to the samples,
 *  and the mixture distribution is updated using expectation-maximization to match the samples.
 */
class AdaptiveImportanceSampler : public BaseSampler {
public:

    /**
     *  @brief Construct a new sampler
     *
     *  @param[in]  rng         Random number generator to use to generate samples.
     *  @param[in]  proposal    Initial distribution to draw from.  It will be modified on each
     *                          each iteration, so users should copy before passing if they
     *                          need to protect the original from modifification.
     *  @param[in]  prior       Bayesian prior used to determine the marginalized posterior
     *                          when updating the proposal distribution to match it.
     *  @param[in]  center      Center position of source.
     *  @param[in]  ctrls       Vector of control objects that define the iterations.
     *  @param[in]  doSaveIterations   Whether to save intermediate SampleSets and associated
     *                                 proposal distributions.
     */
    AdaptiveImportanceSampler(
        PTR(afw::math::Random) rng,
        PTR(MixtureBase) proposal,
        PTR(Prior const) prior,
        afw::geom::Point2D const & center,
        std::map<int,ImportanceSamplerControl> const & ctrls,
        bool doSaveIterations=false
    );

    /**
     *  @brief Generate and evaluate samples using adaptive importance sampling
     */
    virtual SampleSet run(Objective const & objective) const;

    typedef std::map< int, std::vector<PTR(SampleSet)> > IterationMap;

    /**
     *  @brief Return the SampleSet corresponding to the given iteration and repeat numbers.
     *
     *  Only valid if saveIterations=true was passed on initialization.
     *
     *  The proposal distribution attached to the SampleSet will be a snapshot from before to any
     *  E-M updates in that iteration, so the proposal distribution reflects the state of the
     *  proposal the samples were actually drawn from.  This is different from the final SampleSet
     *  returned by the run() method, whose attached proposal distribution does reflect any final
     *  update steps.
     */
    IterationMap const & getIterations() const { return _iterations; }

private:
    bool _doSaveIterations;
    PTR(afw::math::Random)  _rng;
    PTR(MixtureBase) _proposal;
    PTR(Prior const) _prior;
    afw::geom::Point2D _center;
    std::map<int,ImportanceSamplerControl> _ctrls;
    mutable IterationMap _iterations;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_AdaptiveImportanceSampler_h_INCLUDED
