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

#include "ndarray/eigen.h"

#include "lsst/utils/ieee.h"
#include "lsst/afw/table/BaseRecord.h"
#include "lsst/afw/table/Catalog.h"
#include "lsst/meas/multifit/AdaptiveImportanceSampler.h"

namespace lsst { namespace meas { namespace multifit {

namespace {

// Given a sample catalog with log unnormalized weights, transform to normalized weights
Scalar computeRobustWeights(afw::table::BaseCatalog & samples, afw::table::Key<Scalar> const & weightKey) {
    // Sort the sample by weight so we can accumulate robustly.
    samples.sort(weightKey);
    // We now compute z, the arithmetic mean of ln(p_i/q_i).
    Scalar z = 0.0;
    for (afw::table::BaseCatalog::iterator i = samples.begin(); i != samples.end(); ++i) {
        z += i->get(weightKey);
    }
    z /= samples.size();
    // We now subtract z from w_i and exponentiate, accumulating the sums.
    // This makes w_i = e^{-z} p_i/q_i, which is proportional to the
    // desired p_i/q_i.
    Scalar wSum = 0.0;
    for (afw::table::BaseCatalog::iterator i = samples.begin(); i != samples.end(); ++i) {
        wSum += (*i)[weightKey] = std::exp(i->get(weightKey) - z);
        if (!utils::isfinite(i->get(weightKey))) {
            throw LSST_EXCEPT(
                pex::exceptions::LogicErrorException,
                (boost::format(
                    "weight = %g not finite in samples before normalization")
                    % i->get(weightKey)).str()
            );
        }
    }
    if (wSum <= 0.0) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            (boost::format("wSum = %g not positive in samples before normalization") % wSum).str()
        );
    }
    // finally, we normalize w_i...
    for (afw::table::BaseCatalog::iterator i = samples.begin(); i != samples.end(); ++i) {
        (*i)[weightKey] /= wSum;
        if (!utils::isfinite(i->get(weightKey))) {
            throw LSST_EXCEPT(
                pex::exceptions::LogicErrorException,
                (boost::format(
                    "sample weight %g not finite after normalization")
                    % i->get(weightKey)).str()
            );
        }
    }
    // ..and return the log of wSum, corrected for the z term we took out earlier,
    // and including the r/2 term we've ignored all along.
    return - z - std::log(wSum / samples.size());
}

} // anonymous

AdaptiveImportanceSampler::AdaptiveImportanceSampler(
    afw::table::Schema & sampleSchema,
    PTR(afw::math::Random) rng,
    std::map<int,ImportanceSamplerControl> const & ctrls,
    bool doSaveIterations
) :
    _doSaveIterations(doSaveIterations),
    _rng(rng),
    _ctrls(ctrls),
    _weightKey(sampleSchema.addField<Scalar>("weight", "normalized Monte Carlo weight")),
    _objectiveKey(sampleSchema.addField<Scalar>(
                      "objective", "value of the objective function (usually -log posterior)")),
    _proposalKey(sampleSchema.addField<Scalar>(
                     "proposal", "-log value of the proposal function")),
    _parametersKey(sampleSchema["parameters"])
{
    if (_doSaveIterations) {
        _iterCtrlKey = sampleSchema.addField<int>(
            "iter.ctrl", "iteration major ID; corresponds to control map key"
        );
        _iterRepeatKey = sampleSchema.addField<int>(
            "iter.repeat", "iteration minor ID; corresponds to repeats within a single control object"
        );
    }
}

void AdaptiveImportanceSampler::run(
    SamplerObjective const & objective,
    PTR(MixtureBase) proposal,
    afw::table::BaseCatalog & samples
) const {
    double perplexity = 0.0;
    int parameterDim = objective.getParameterDim();
    for (std::map<int,ImportanceSamplerControl>::const_iterator i = _ctrls.begin(); i != _ctrls.end(); ++i) {
        ImportanceSamplerControl const & ctrl = i->second;
        int nRepeat = 0;
        while (nRepeat <= ctrl.maxRepeat && perplexity < ctrl.targetPerplexity) {
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
                PTR(afw::table::BaseRecord) record = samples.addNew();
                double objectiveValue = objective(parameters[k], *record);
                if (utils::isfinite(objectiveValue)) {
                    subSamples.push_back(record);
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
            computeRobustWeights(subSamples, _weightKey);
            perplexity = computeNormalizedPerplexity(subSamples);
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
        h -= s->get(_weightKey) * std::log(s->get(_weightKey));
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

}}} // namespace lsst::meas::multifit
