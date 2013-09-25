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
#include "lsst/meas/multifit/AdaptiveImportanceSampler.h"

namespace lsst { namespace meas { namespace multifit {

AdaptiveImportanceSampler::AdaptiveImportanceSampler(
    PTR(afw::math::Random) rng,
    PTR(MixtureBase) proposal,
    PTR(Prior const) prior,
    afw::geom::Point2D const & center,
    std::map<int,ImportanceSamplerControl> const & ctrls,
    bool doSaveIterations
) :
    _doSaveIterations(doSaveIterations),
    _rng(rng),
    _proposal(proposal),
    _prior(prior),
    _center(center),
    _ctrls(ctrls)
{}

SampleSet AdaptiveImportanceSampler::run(Likelihood const & likelihood) const {
    // TODO: check that nonlinearDim agrees with ParameterDefinition
    int const nonlinearDim = _proposal->getDimension();
    int const linearDim = likelihood.getLinearDim();
    PTR(ParameterDefinition const) parameterDef
        = ParameterDefinition::makeEllipseCoreDefinition("SeparableConformalShearLogTraceRadius", _center);
    SampleSetKeys keys(nonlinearDim, linearDim);
    SampleSet samples(parameterDef, linearDim);
    samples.setDataSquaredNorm(likelihood.getDataSquaredNorm());
    double perplexity = 0.0;
    std::vector<PTR(SampleSet)> * iterationVector = 0;
    for (std::map<int,ImportanceSamplerControl>::const_iterator i = _ctrls.begin(); i != _ctrls.end(); ++i) {
        ImportanceSamplerControl const & ctrl = i->second;
        int nRepeat = 0;
        if (_doSaveIterations) {
            IterationMap::iterator j = _iterations.insert(
                _iterations.end(),
                std::make_pair(i->first, std::vector<PTR(SampleSet)>())
            );
            iterationVector = &j->second;
        }
        while (nRepeat <= ctrl.maxRepeat && perplexity < ctrl.targetPerplexity) {
            ++nRepeat;
            samples.clear();
            ndarray::Array<samples::Scalar,2,2> parameters = ndarray::allocate(ctrl.nSamples, nonlinearDim);
            _proposal->draw(*_rng, parameters);
            ndarray::Array<samples::Scalar,1,1> probability = ndarray::allocate(ctrl.nSamples);
            _proposal->evaluate(parameters, probability);
            for (int k = 0; k < ctrl.nSamples; ++k) {
                LogGaussian joint = likelihood.evaluate(parameterDef->makeEllipse(parameters[k]));
                samples.add(
                    joint,
                    -std::log(probability[k]),
                    parameters[k].asEigen()
                );
            }
            if (_doSaveIterations || ctrl.nUpdateSteps > 0 || nRepeat <= ctrl.maxRepeat) {
                samples.applyPrior(_prior);
                perplexity = samples.computeNormalizedPerplexity();
            }
            if (_doSaveIterations) {
                PTR(SampleSet) snapshot(new SampleSet(samples));
                iterationVector->push_back(snapshot);
                snapshot->setProposal(_proposal->clone());
            }
            if (ctrl.nUpdateSteps > 0) {
                afw::table::BaseCatalog cat = samples.getCatalog();
                for (std::size_t k = 0; k < cat.size(); ++k) {
                    parameters[k] = cat[k].get(keys.parameters);
                    probability[k] = cat[k].get(keys.weight);
                }
                for (int j = 0; j < ctrl.nUpdateSteps; ++j) {
                    _proposal->updateEM(
                        parameters[ndarray::view(0, cat.size())],
                        probability[ndarray::view(0, cat.size())],
                        ctrl.tau1, ctrl.tau2
                    );
                }
            }
        }
    }
    samples.setProposal(_proposal->clone());
    return samples;
}

}}} // namespace lsst::meas::multifit
