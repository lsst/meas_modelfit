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

#include "lsst/meas/modelfit/ModelFitRecord.h"
#include "lsst/meas/modelfit/Sampling.h"

namespace lsst { namespace meas { namespace modelfit {

ndarray::Array<Scalar,1,1> SamplingInterpreter::computeParameterQuantiles(
    ModelFitRecord const & record,
    ndarray::Array<Scalar const,1,1> const & fractions,
    int index
) const {
    return computeSampleQuantiles(record.getSamples(), fractions, _parameterKey[index]);
}

ndarray::Array<Scalar,1,1> SamplingInterpreter::computeNonlinearQuantiles(
    ModelFitRecord const & record,
    ndarray::Array<Scalar const,1,1> const & fractions,
    int index
) const {
    return computeSampleQuantiles(record.getSamples(), fractions, _nonlinearKey[index]);
}

ndarray::Array<Scalar,1,1> SamplingInterpreter::computeParameterMean(
    ModelFitRecord const & record
) const {
    return computeSampleMean(record.getSamples(), _parameterKey);
}

ndarray::Array<Scalar,1,1> SamplingInterpreter::computeNonlinearMean(
    ModelFitRecord const & record
) const {
    return computeSampleMean(record.getSamples(), _nonlinearKey);
}

ndarray::Array<Scalar,2,2> SamplingInterpreter::computeParameterCovariance(
    ModelFitRecord const & record,
    ndarray::Array<Scalar const,1,1> const & mean
) const {
    return computeSampleCovariance(record.getSamples(), mean, _parameterKey);
}

ndarray::Array<Scalar,2,2> SamplingInterpreter::computeNonlinearCovariance(
    ModelFitRecord const & record,
    ndarray::Array<Scalar const,1,1> const & mean
) const {
    return computeSampleCovariance(record.getSamples(), mean, _nonlinearKey);
}

SamplingInterpreter::SamplingInterpreter(
    afw::table::Schema & sampleSchema,
    Model::NameVector const & parameterNames,
    PTR(Model) model,
    PTR(Prior) prior
) : Interpreter(parameterNames, model, prior),
    _weightKey(
        sampleSchema.addField(
            afw::table::Field<Scalar>("weight", "normalized Monte Carlo weight"),
            true // doReplace
        )
    )
{}

ndarray::Array<Scalar,1,1> SamplingInterpreter::computeSampleQuantiles(
    afw::table::BaseCatalog const & samples,
    ndarray::Array<Scalar const,1,1> const & fractions,
    ScalarKey const & targetKey
) const {
    // in all std::pairs below, first == parameter value, second == (possibly cumulative) weight
    ndarray::Array<Scalar,1,1> result = ndarray::allocate(fractions.getSize<0>());
    result.deep() = 0.0;
    if (!samples.empty()) {
        // use map to sort by parameter value, merge entries with exactly equal parameter values
        std::map<Scalar,Scalar> map;
        for (afw::table::BaseCatalog::const_iterator i = samples.begin(); i != samples.end(); ++i) {
            std::pair<std::map<Scalar,Scalar>::iterator,bool> r = map.insert(
                std::pair<Scalar,Scalar>(i->get(targetKey), i->get(_weightKey))
            );
            if (!r.second) r.first->second += i->get(_weightKey);
        }
        Scalar cumulative = 0.0;
        std::size_t iFraction = 0;
        std::map<Scalar,Scalar>::const_iterator current = map.begin(), end = map.end();
        std::pair<Scalar,Scalar> last(current->first, 0.0);
        for (; current != end; ++current) {
            cumulative += current->second;
            // see if we exceeded one of the desired fractions; if so, we
            // linearly interpolate the exact point at which we would have
            // met that fraction
            while (cumulative >= fractions[iFraction]) {
                Scalar delta = current->first - last.first;
                Scalar w = (fractions[iFraction] - last.second) / current->second;
                result[iFraction] = current->first + w * delta;
                if (result[iFraction] > current->first) {  // can happen due to round-off error
                    result[iFraction] = current->first;
                }
                ++iFraction;
                if (iFraction == fractions.size()) break;
            }
            if (iFraction == fractions.size()) break;
            last.first = current->first;
            last.second = cumulative;
        }
        while (iFraction < fractions.size()) {
            result[iFraction] = last.first;
            ++iFraction;
        }
    }
    return result;
}

ndarray::Array<Scalar,1,1> SamplingInterpreter::computeSampleMean(
    afw::table::BaseCatalog const & samples,
    ArrayKey const & targetKey
) const {
    ndarray::Array<Scalar,1,1> result = ndarray::allocate(targetKey.getSize());
    result.deep() = 0.0;
    Scalar wSum = 0.0;
    for (afw::table::BaseCatalog::const_iterator i = samples.begin(); i != samples.end(); ++i) {
        wSum += i->get(_weightKey);
        result.asEigen() += i->get(targetKey).asEigen() * i->get(_weightKey);
    }
    result.asEigen() /= wSum;
    return result;
}

ndarray::Array<Scalar,2,2> SamplingInterpreter::computeSampleCovariance(
    afw::table::BaseCatalog const & samples,
    ndarray::Array<Scalar const,1,1> const & mean,
    ArrayKey const & targetKey
) const {
    ndarray::Array<Scalar,2,2> result = ndarray::allocate(targetKey.getSize(), targetKey.getSize());
    result.deep() = 0.0;
    Scalar wSum = 0.0;
    Vector workspace;
    for (afw::table::BaseCatalog::const_iterator i = samples.begin(); i != samples.end(); ++i) {
        workspace = i->get(targetKey).asEigen() - mean.asEigen();
        wSum += i->get(_weightKey);
        result.asEigen().selfadjointView<Eigen::Lower>().rankUpdate(workspace, i->get(_weightKey));
    }
    result.asEigen() = result.asEigen().selfadjointView<Eigen::Lower>();
    result.asEigen() /= wSum;
    return result;
}

SamplingObjective::SamplingObjective(
    PTR(SamplingInterpreter) interpreter,
    PTR(Likelihood) likelihood
) :
    _interpreter(interpreter),
    _likelihood(likelihood),
    _modelMatrix(ndarray::allocate(likelihood->getDataDim(), likelihood->getAmplitudeDim()))
{}

}}} // namespace lsst::meas::modelfit
