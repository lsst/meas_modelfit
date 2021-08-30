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

#ifndef LSST_MEAS_MODELFIT_AdaptiveImportanceSampler_h_INCLUDED
#define LSST_MEAS_MODELFIT_AdaptiveImportanceSampler_h_INCLUDED

#include <map>

#include "lsst/pex/config.h"
#include "lsst/afw/table/Schema.h"
#include "lsst/meas/modelfit/Sampler.h"
#include "lsst/meas/modelfit/Mixture.h"

namespace lsst { namespace meas { namespace modelfit {

/**
 *  @brief Control object for one iteration of adaptive importance sampling
 *
 *  @sa AdaptiveImportanceSampler, AdaptiveImportanceSamplerTask
 */
class ImportanceSamplerControl {
public:
    LSST_CONTROL_FIELD(nSamples, int, "Number of Monte Carlo samples to draw");
    LSST_CONTROL_FIELD(nUpdateSteps, int, "Number of Expectation-Maximization update iterations");
    LSST_CONTROL_FIELD(tau1, double, "Damping parameter for E-M update (see Mixture::updateEM)");
    LSST_CONTROL_FIELD(tau2, double, "Damping parameter for E-M update (see Mixture::updateEM)");
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
class AdaptiveImportanceSampler : public Sampler {
public:

    /**
     *  @brief Construct a new sampler
     *
     *  @param[in,out] sampleSchema   Schema for the catalog of samples filled by the Sampler;
     *                                will be modified to include sampler-specific fields.
     *  @param[in]  rng               Random number generator to use to generate samples.
     *  @param[in]  ctrls             Vector of control objects that define the iterations.
     *  @param[in]  doSaveIterations  Whether to save intermediate SampleSets and associated
     *                                proposal distributions.
     */
    AdaptiveImportanceSampler(
        afw::table::Schema & sampleSchema,
        std::shared_ptr<afw::math::Random> rng,
        std::map<int,ImportanceSamplerControl> const & ctrls,
        bool doSaveIterations=false
    );

    void run(
        SamplingObjective const & objective,
        std::shared_ptr<Mixture> proposal,
        afw::table::BaseCatalog & samples
    ) const override;

    double computeNormalizedPerplexity(afw::table::BaseCatalog const & samples) const;

    double computeEffectiveSampleSizeFraction(afw::table::BaseCatalog const & samples) const;

private:
    bool _doSaveIterations;
    std::shared_ptr<afw::math::Random>  _rng;
    std::map<int,ImportanceSamplerControl> _ctrls;
    afw::table::Key<Scalar> _weightKey;
    afw::table::Key<Scalar> _objectiveKey;
    afw::table::Key<Scalar> _proposalKey;
    afw::table::Key< afw::table::Array<Scalar> > _parametersKey;
    afw::table::Key<int> _iterCtrlKey;
    afw::table::Key<int> _iterRepeatKey;
};

}}} // namespace lsst::meas::modelfit

#endif // !LSST_MEAS_MODELFIT_AdaptiveImportanceSampler_h_INCLUDED
