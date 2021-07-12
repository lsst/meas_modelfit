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
#include <cmath>

#include "ndarray/eigen.h"

#include "lsst/log/Log.h"
#include "lsst/afw/table/BaseRecord.h"
#include "lsst/afw/table/Catalog.h"
#include "lsst/meas/modelfit/AdaptiveImportanceSampler.h"

namespace lsst { namespace meas { namespace modelfit {

namespace {

// Given a sample catalog with log unnormalized weights, transform to normalized weights
Scalar computeRobustWeights(afw::table::BaseCatalog & samples, afw::table::Key<Scalar> const & weightKey) {
    LOG_LOGGER trace4Logger = LOG_GET("TRACE4.meas.modelfit.AdaptiveImportanceSampler");
    static Scalar const CLIP_THRESHOLD = 100; // clip samples with weight < e^{-CLIP_THRESHOLD} * wMax
    LOGL_DEBUG(trace4Logger, "Starting computeRobustWeights with %d samples", int(samples.size()));
    // Sort the sample by weight so we can accumulate robustly.
    samples.sort(weightKey);
    Scalar uMax = samples.back().get(weightKey);
    Scalar uClip = uMax - CLIP_THRESHOLD;
    afw::table::BaseCatalog::iterator const iClip = samples.lower_bound(uClip, weightKey);
    LOGL_DEBUG(trace4Logger, "uMax=%g, uClip=%g, iClip at offset %d", uMax, uClip, int(iClip - samples.begin()));
    for (afw::table::BaseCatalog::iterator i = samples.begin(); i != iClip; ++i) {
        i->set(weightKey, 0.0);
    }
    Scalar wSum = 0.0;
    for (afw::table::BaseCatalog::iterator i = iClip; i != samples.end(); ++i) {
        Scalar w = std::exp(i->get(weightKey) - uMax);
        i->set(weightKey, w);
        wSum += w;
    }
    LOGL_DEBUG(trace4Logger, "Uncorrected wSum=%g", wSum);
    for (afw::table::BaseCatalog::iterator i = iClip; i != samples.end(); ++i) {
        (*i)[weightKey] /= wSum;
    }
    return - uMax - std::log(wSum / samples.size());
}

} // anonymous

AdaptiveImportanceSampler::AdaptiveImportanceSampler(
    afw::table::Schema & sampleSchema,
    std::shared_ptr<afw::math::Random> rng,
    std::map<int,ImportanceSamplerControl> const & ctrls,
    bool doSaveIterations
) :
    _doSaveIterations(doSaveIterations),
    _rng(rng),
    _ctrls(ctrls),
    _weightKey(sampleSchema["weight"]),
    _objectiveKey(
        sampleSchema.addField(
            afw::table::Field<Scalar>(
                "objective", "value of the objective function (usually -log posterior)"
            ),
            true // doReplace
        )
    ),
    _proposalKey(
        sampleSchema.addField(
            afw::table::Field<Scalar>("proposal", "-log value of the proposal function"),
            true // doReplace
        )
    ),
    _parametersKey(sampleSchema["parameters"])
{
    if (_doSaveIterations) {
        _iterCtrlKey = sampleSchema.addField(
            afw::table::Field<int>(
                "iter.ctrl", "iteration major ID; corresponds to control map key"
            ),
            true // doReplace
        );
        _iterRepeatKey = sampleSchema.addField(
            afw::table::Field<int>(
                "iter.repeat", "iteration minor ID; corresponds to repeats within a single control object"
            ),
            true // doReplace
        );
    }
}

void AdaptiveImportanceSampler::run(
    SamplingObjective const & objective,
    std::shared_ptr<Mixture> proposal,
    afw::table::BaseCatalog & samples
) const {
    LOG_LOGGER trace3Logger = LOG_GET("TRACE3.meas.modelfit.AdaptiveImportanceSampler");
    double perplexity = 0.0;
    int parameterDim = objective.getParameterDim();
    for (std::map<int,ImportanceSamplerControl>::const_iterator i = _ctrls.begin(); i != _ctrls.end(); ++i) {
        ImportanceSamplerControl const & ctrl = i->second;
        int nRepeat = 0;
        while (nRepeat <= ctrl.maxRepeat && perplexity < ctrl.targetPerplexity) {
            LOGL_DEBUG(trace3Logger,
                "Starting repeat %d with nSamples=%d, nUpdateSteps=%d, targetPerplexity=%g",
                nRepeat, ctrl.nSamples, ctrl.nUpdateSteps, ctrl.targetPerplexity
            );
            ++nRepeat;
            if (!_doSaveIterations) {
                samples.clear();
            }
            afw::table::BaseCatalog subSamples(samples.getTable());
            ndarray::Array<Scalar,2,2> parameters = ndarray::allocate(ctrl.nSamples, parameterDim);
            proposal->draw(*_rng, parameters);
            ndarray::Array<Scalar,1,1> probability = ndarray::allocate(ctrl.nSamples);
            proposal->evaluate(parameters, probability);
            for (int k = 0; k < ctrl.nSamples; ++k) {
                std::shared_ptr<afw::table::BaseRecord> record = samples.addNew();
                double objectiveValue = objective(parameters[k], *record);
                if (std::isfinite(objectiveValue)) {
                    subSamples.push_back(record);
                    record->set(_parametersKey, parameters[k]);
                    record->set(_objectiveKey, objectiveValue);
                    record->set(_proposalKey, -std::log(probability[k]));
                    if (_doSaveIterations) {
                        record->set(_iterCtrlKey, i->first);
                        record->set(_iterRepeatKey, nRepeat-1);
                    }
                    // for numerical reasons, in the first pass, we set w_i = ln(p_i/q_i);
                    // note that proposal[i] == -ln(q_i) and objective[i] == -ln(p_i)
                    record->set(_weightKey, record->get(_proposalKey) - record->get(_objectiveKey));
                } else {
                    samples.pop_back();
                }
            }
            if (samples.empty()) {
                throw LSST_EXCEPT(
                    pex::exceptions::LogicError,
                    "No finite objective values in entire sample set"
                );
            }
            computeRobustWeights(subSamples, _weightKey);
            perplexity = computeNormalizedPerplexity(subSamples);
            if (!std::isfinite(perplexity)) {
                throw LSST_EXCEPT(
                    pex::exceptions::LogicError,
                    "Normalized perplexity is non-finite."
                );
            }
            LOGL_DEBUG(trace3Logger,
                "Normalized perplexity is %g; target is %g",
                perplexity, ctrl.targetPerplexity
            );
            if (ctrl.nUpdateSteps > 0) {
                for (std::size_t k = 0; k < subSamples.size(); ++k) {
                    parameters[k] = subSamples[k].get(_parametersKey);
                    probability[k] = subSamples[k].get(_weightKey);
                }
                for (int j = 0; j < ctrl.nUpdateSteps; ++j) {
                    proposal->updateEM(
                        parameters[ndarray::view(0, subSamples.size())],
                        probability[ndarray::view(0, subSamples.size())],
                        ctrl.tau1, ctrl.tau2
                    );
                }
            }
        }
    }
}

double AdaptiveImportanceSampler::computeNormalizedPerplexity(
    afw::table::BaseCatalog const & samples
) const {
    double h = 0.0;
    for (afw::table::BaseCatalog::const_iterator s = samples.begin(); s != samples.end(); ++s) {
        if (s->get(_weightKey) > 0.0) {
            h -= s->get(_weightKey) * std::log(s->get(_weightKey));
        }
    }
    return std::exp(h) / samples.size();
}

double AdaptiveImportanceSampler::computeEffectiveSampleSizeFraction(
    afw::table::BaseCatalog const & samples
) const {
    double t = 0.0;
    for (afw::table::BaseCatalog::const_iterator s = samples.begin(); s != samples.end(); ++s) {
        t += s->get(_weightKey) * s->get(_weightKey);
    }
    return 1.0 / (t * samples.size());
}

}}} // namespace lsst::meas::modelfit
