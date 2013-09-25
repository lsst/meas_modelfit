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

#include <limits>

#include "ndarray/eigen.h"

#include "lsst/utils/ieee.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/table/BaseColumnView.h"
#include "lsst/afw/table/io/CatalogVector.h"
#include "lsst/afw/table/io/InputArchive.h"
#include "lsst/afw/table/io/OutputArchive.h"
#include "lsst/meas/multifit/SampleSet.h"
#include "lsst/meas/multifit/ExpectationFunctor.h"
#include "lsst/meas/multifit/priors.h"

namespace tbl = lsst::afw::table;

namespace lsst { namespace meas { namespace multifit {

SampleSetKeys::SampleSetKeys(int nonlinearDim, int linearDim) :
    schema(),
    jointGrad(schema.addField<samples::ArrayTag>(
                "joint.grad", "amplitude log-likelihood gradient at amplitude=0", linearDim)),
    jointFisher(schema.addField<samples::ArrayTag>("joint.fisher", "amplitude Fisher matrix",
                                                   linearDim*linearDim)),
    marginal(schema.addField<samples::Scalar>(
                 "marginal", "negative log marginal posterior value at sample point")),
    proposal(schema.addField<samples::Scalar>(
                 "proposal", "negative log density of the distribution used to draw samples")),
    weight(schema.addField<samples::Scalar>("weight", "normalized Monte Carlo weight")),
    parameters(schema.addField<samples::ArrayTag>(
                   "parameters", "nonlinear parameters at this point", nonlinearDim))
{}

SampleSetKeys::SampleSetKeys(tbl::Schema const & schema_) :
    schema(schema_),
    jointGrad(schema["joint.grad"]),
    jointFisher(schema["joint.fisher"]),
    marginal(schema["marginal"]),
    proposal(schema["proposal"]),
    weight(schema["weight"]),
    parameters(schema["parameters"])
{}

LogGaussian SampleSetKeys::getJoint(tbl::BaseRecord const & record) const {
    LogGaussian joint(getLinearDim());
    joint.grad = samples::VectorCMap(record.getElement(jointGrad), getLinearDim());
    joint.fisher = samples::MatrixCMap(record.getElement(jointFisher), getLinearDim(), getLinearDim());
    return joint;
}

void SampleSetKeys::setJoint(tbl::BaseRecord & record, LogGaussian const & joint) {
    samples::VectorMap(record.getElement(jointGrad), getLinearDim()) = joint.grad;
    samples::MatrixMap(record.getElement(jointFisher), getLinearDim(), getLinearDim()) = joint.fisher;
    if (*record.getElement(jointFisher) != joint.fisher(0,0)) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "*record.getElement(jointFisher) != joint.fisher(0,0) in SampleSetKeys::setJoint"
        );
    }
}

SampleSet::SampleSet(PTR(ParameterDefinition const) parameterDefinition, int linearDim) :
    _keys(parameterDefinition->getDim(), linearDim),
    _records(_keys.schema),
    _dataSquaredNorm(0.0),
    _parameterDefinition(parameterDefinition),
    _prior()
{}

SampleSet::SampleSet(PTR(ParameterDefinition const) parameterDefinition, tbl::BaseCatalog const & records) :
    _keys(records.getSchema()),
    _records(records),
    _dataSquaredNorm(0.0),
    _parameterDefinition(parameterDefinition),
    _prior()
{}

void SampleSet::setParameterDefinition(PTR(ParameterDefinition const) parameterDefinition) {
    if (*_parameterDefinition != *parameterDefinition) {
        PTR(ParameterConverter const) converter = _parameterDefinition->makeConverterTo(*parameterDefinition);
        for (tbl::BaseCatalog::iterator s = _records.begin(); s != _records.end(); ++s) {
            ndarray::Array<samples::Scalar,1,1> a = (*s)[_keys.parameters];
            converter->apply(a, a);
        }
    }
    _parameterDefinition = parameterDefinition;
}

ndarray::Array<double,1,1> SampleSet::computeDensity(KernelDensityEstimatorControl const & ctrl) const {
    tbl::Key<samples::Scalar> paramKey = _keys.parameters[ctrl.getIndex()];
    ndarray::Array<double,1,1> result = ndarray::allocate(ctrl.getCount());
    result.deep() = 0.0;
    for (tbl::BaseCatalog::const_iterator s = _records.begin(); s != _records.end(); ++s) {
        ctrl.apply(result, s->get(paramKey), s->get(_keys.weight));
    }
    return result;
}

ndarray::Array<double,2,2> SampleSet::computeDensity(
    KernelDensityEstimatorControl const & ctrlX,
    KernelDensityEstimatorControl const & ctrlY
) const {
    tbl::Key<samples::Scalar> paramKeyX = _keys.parameters[ctrlX.getIndex()];
    tbl::Key<samples::Scalar> paramKeyY = _keys.parameters[ctrlY.getIndex()];
    ndarray::Array<double,2,2> result = ndarray::allocate(ctrlY.getCount(), ctrlX.getCount());
    ndarray::Array<double,1,1> tmpY = ndarray::allocate(ctrlY.getCount());
    ndarray::Array<double,1,1> tmpX = ndarray::allocate(ctrlX.getCount());
    result.deep() = 0.0;
    for (tbl::BaseCatalog::const_iterator s = _records.begin(); s != _records.end(); ++s) {
        tmpY.deep() = 0.0;
        tmpX.deep() = 0.0;
        ctrlY.apply(tmpY, s->get(paramKeyY), 1.0);
        ctrlX.apply(tmpX, s->get(paramKeyX), 1.0);
        result.asEigen() += s->get(_keys.weight) * tmpY.asEigen() * tmpX.asEigen().adjoint();
    }
    return result;
}

double SampleSet::computeNormalizedPerplexity() const {
    double h = 0.0;
    for (tbl::BaseCatalog::const_iterator s = _records.begin(); s != _records.end(); ++s) {
        h -= s->get(_keys.weight) * std::log(s->get(_keys.weight));
    }
    return std::exp(h) / _records.size();
}

double SampleSet::computeEffectiveSampleSizeFraction() const {
    double t = 0.0;
    for (tbl::BaseCatalog::const_iterator s = _records.begin(); s != _records.end(); ++s) {
        t += s->get(_keys.weight) * s->get(_keys.weight);
    }
    return 1.0 / (t * _records.size());
}

void SampleSet::add(LogGaussian const & joint, double proposal, samples::Vector const & parameters) {
    if (parameters.size() != getNonlinearDim()) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            (boost::format("Incorrect nonlinear dimension for SamplePoint: expected %d, got %d")
             % getNonlinearDim() % parameters.size()).str()
        );
    }
    if (joint.grad.size() != getLinearDim()) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            (boost::format("Incorrect linear dimension for SamplePoint: expected %d, got %d")
             % getNonlinearDim() % joint.grad.size()).str()
        );
    }
    if (_prior) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "Cannot add new samples after applyPrior() has been called; use dropPrior() to revert"
        );
    }
    PTR(tbl::BaseRecord) newRecord = _records.addNew();
    newRecord->set(_keys.proposal, proposal);
    _keys.setJoint(*newRecord, joint);
    _keys.setParameters(*newRecord, parameters);
}

double SampleSet::applyPrior(PTR(Prior const) prior, double clip) {
    _prior = prior;
    std::size_t origSize = _records.size();
    tbl::BaseCatalog::iterator i = _records.begin();
    ndarray::Array<double,1,1> priorParameters;
    PTR(ParameterConverter const) converter;
    if (prior->getParameterDefinition() != getParameterDefinition()) {
        converter = getParameterDefinition()->makeConverterTo(*prior->getParameterDefinition());
        priorParameters = ndarray::allocate(getNonlinearDim());
    }
    double logJacobian = 0.0;
    for (; i != _records.end(); ++i) {
        if (converter) {
            converter->apply(i->get(_keys.parameters), priorParameters);
            logJacobian = std::log(converter->computeJacobian(i->get(_keys.parameters)));
        } else {
            priorParameters = (*i)[_keys.parameters];
        }
        i->set(_keys.marginal, prior->apply(_keys.getJoint(*i), priorParameters.asEigen()) + logJacobian);
        if (!utils::isfinite(logJacobian)) {
            // Sometimes we get NaNs in Jacobian due to 1/r where r==0; these are actually infinities,
            // and should translate into zero-probability points.
            i->set(_keys.marginal, std::numeric_limits<samples::Scalar>::infinity());
        }
        // for numerical reasons, in the first pass, we set w_i = ln(m_i/q_i);
        // note that i->proposal == -ln(q_i) and i->marginal == -ln(m_i)
        i->set(_keys.weight, i->get(_keys.proposal) - i->get(_keys.marginal));
        if (utils::isnan(i->get(_keys.weight))) {
            throw LSST_EXCEPT(
                pex::exceptions::RuntimeErrorException,
                (boost::format("NaN encountered in weights: marginal=%f, proposal=%f")
                 % i->get(_keys.marginal) % i->get(_keys.proposal)).str()
            );
        }
    }
    // sort by ascending probability, so when we accumulate, we add small numbers together
    // before adding them to large numbers
    _records.sort(_keys.weight);
    // now we clip weights that are much smaller than the maximum weight, as we know they're
    // insignificant in the full sum, and they can cause numerical issues later on
    clip = std::max(clip, std::numeric_limits<double>::min());
    double const tau = std::log(clip) + _records.back().get(_keys.weight);
    i = _records.begin();
    while (i != _records.end()) {
        if (i->get(_keys.weight) >= tau) break;
        ++i;
    }
    _records.erase(_records.begin(), i);
    if (_records.empty()) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            "No unclipped samples after applyPrior"
        );
    }
    // we now compute z, the arithmetic mean of ln(m_i/q_i)
    double z = 0.0;
    for (i = _records.begin(); i != _records.end(); ++i) {
        z += i->get(_keys.weight);
    }
    z /= _records.size();
    // we now subtract z from w_i and exponentiate, accumulating the sums
    // this now makes w_i = e^{-z} m_i/q_i, which is proportional to the
    // desired m_i/q_i.
    double wSum = 0.0;
    for (i = _records.begin(); i != _records.end(); ++i) {
        wSum += (*i)[_keys.weight] = std::exp(i->get(_keys.weight) - z);
        if (!utils::isfinite(i->get(_keys.weight))) {
            throw LSST_EXCEPT(
                pex::exceptions::LogicErrorException,
                (boost::format(
                    "i->get(_keys.weight) = %d not finite in SampleSet::applyPrior before normalization")
                    % i->get(_keys.weight)).str()
            );
        }
    }
    if (wSum <= 0.0) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            (boost::format("wSum = %d not positive in SampleSet::applyPrior") % wSum).str()
        );
    }
    // finally, we normalize w_i...
    for (i = _records.begin(); i != _records.end(); ++i) {
        (*i)[_keys.weight] /= wSum;
        if (!utils::isfinite(i->get(_keys.weight))) {
            throw LSST_EXCEPT(
                pex::exceptions::LogicErrorException,
                (boost::format(
                    "i->get(_keys.weight) = %d not finite in SampleSet::applyPrior after normalization")
                    % i->get(_keys.weight)).str()
            );
        }
    }
    // ..and return the log of wSum, corrected for the z term we took out earlier,
    // and including the r/2 term we've ignored all along.
    return 0.5*_dataSquaredNorm - z - std::log(wSum / origSize);
}

void SampleSet::dropPrior() {
    for (tbl::BaseCatalog::iterator i = _records.begin(); i != _records.end(); ++i) {
        (*i)[_keys.weight] = (*i)[_keys.marginal] = std::numeric_limits<double>::quiet_NaN();
    }
    _prior.reset();
}

void SampleSet::clear() {
    _records.clear();
    _prior.reset();
}

samples::Vector SampleSet::computeQuantiles(samples::Vector const & fraction, int parameterIndex) const {
    samples::ScalarKey paramKey = _keys.parameters[parameterIndex];
    // in all std::pairs below, first == parameter value, second == (possibly cumulative) weight
    samples::Vector result = samples::Vector::Zero(fraction.size());
    if (!_records.empty()) {
        // use map to sort by parameter value, merge entries with exactly equal parameter values
        std::map<double,double> map;
        for (tbl::BaseCatalog::const_iterator i = _records.begin(); i != _records.end(); ++i) {
            std::pair<std::map<double,double>::iterator,bool> r = map.insert(
                std::pair<double,double>(i->get(paramKey), i->get(_keys.weight))
            );
            if (!r.second) r.first->second += i->get(_keys.weight);
        }
        double cumulative = 0.0;
        int iFraction = 0;
        std::map<double,double>::const_iterator current = map.begin(), end = map.end();
        std::pair<double,double> last(current->first, 0.0);
        for (; current != end; ++current) {
            cumulative += current->second;
            while (cumulative >= fraction[iFraction]) {
                double delta = current->first - last.first;
                double w = (fraction[iFraction] - last.second) / current->second;
                result[iFraction] = current->first + w * delta;
                if (result[iFraction] > current->first) {  // can happen due to round-off error
                    result[iFraction] = current->first;
                }
                ++iFraction;
                if (iFraction == fraction.size()) break;
            }
            if (iFraction == fraction.size()) break;
            last.first = current->first;
            last.second = cumulative;
        }
        while (iFraction < fraction.size()) {
            result[iFraction] = last.first;
            ++iFraction;
        }
    }
    return result;
}

samples::Matrix SampleSet::computeQuantiles(samples::Vector const & fractions) const {
    samples::Matrix result(fractions.size(), getNonlinearDim());
    for (int i = 0; i < getNonlinearDim(); ++i) {
        result.col(i) = computeQuantiles(fractions, i);
    }
    return result;
}

namespace {

class MeanExpectationFunctor : public ExpectationFunctor {
public:
    explicit MeanExpectationFunctor(int dim_) : ExpectationFunctor(dim_) {}

    virtual samples::Vector operator()(samples::VectorCMap const & parameters) const {
        return parameters;
    }
};

class CovarianceExpectationFunctor : public ExpectationFunctor {
public:
    explicit CovarianceExpectationFunctor(samples::Vector const & mean) :
        ExpectationFunctor(mean.size() * (mean.size() - 1)), _mean(mean)
    {}

    virtual samples::Vector operator()(samples::VectorCMap const & parameters) const {
        samples::Vector r = samples::Vector::Zero(outputDim);
        int n = 0;
        for (int i = 0; i < _mean.size(); ++i) {
            for (int j = 0; j <= i; ++j) {
                r[n] = (parameters[i] - _mean[i]) * (parameters[j] - _mean[j]);
                ++n;
            }
        }
        return r;
    }

    samples::Matrix expand(samples::Vector const & x) const {
        samples::Matrix r = samples::Matrix::Zero(_mean.size(), _mean.size());
        int n = 0;
        for (int i = 0; i < _mean.size(); ++i) {
            for (int j = 0; j <= i; ++j) {
                r(j, i) = r(i, j) = x[n];
                ++n;
            }
        }
        return r;
    }

private:
    samples::Vector const & _mean;
};

} // anonymous

samples::Vector SampleSet::computeExpectation(
    ExpectationFunctor const & functor,
    samples::Matrix * mcCov
) const {
    if (!_prior) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "Cannot compute expectation values without attaching a prior"
        );
    }
    samples::Vector r = samples::Vector::Zero(functor.outputDim);
    for (tbl::BaseCatalog::const_iterator i = _records.begin(); i != _records.end(); ++i) {
        r += i->get(_keys.weight) * functor(_keys.getParameters(*i));
    }
    if (mcCov) {
        mcCov->setZero();
        double w2 = 0.0;
        for (tbl::BaseCatalog::const_iterator i = _records.begin(); i != _records.end(); ++i) {
            samples::Vector delta = functor(_keys.getParameters(*i));
            delta -= r;
            mcCov->selfadjointView<Eigen::Lower>().rankUpdate(delta, i->get(_keys.weight));
            w2 += i->get(_keys.weight) * i->get(_keys.weight);
        }
        *mcCov = mcCov->selfadjointView<Eigen::Lower>();
        *mcCov /= (1.0 - w2);
    }
    return r;
}

samples::Vector SampleSet::computeMean(samples::Matrix * mcCov) const {
    MeanExpectationFunctor f(getNonlinearDim());
    return computeExpectation(f, mcCov);
}

samples::Matrix SampleSet::computeCovariance(samples::Vector const & mean) const {
    CovarianceExpectationFunctor f(mean);
    return f.expand(computeExpectation(f));
}

namespace {

// schema and keys for single-record (per SampleSet), dimension-independent catalog,
// used to persist all the SampleSet data members besides the actual catalog of samples
class SampleSetPersistenceKeys : private boost::noncopyable {
public:
    tbl::Schema schema;
    tbl::Key<int> prior;
    tbl::Key<int> proposal;
    tbl::Key<double> dataSquaredNorm;
    tbl::Key<int> parameterDefinition;

    static SampleSetPersistenceKeys const & get() {
        static SampleSetPersistenceKeys const instance;
        return instance;
    }

private:

    SampleSetPersistenceKeys() :
        schema(),
        prior(schema.addField<int>("prior", "archive ID for Bayesian prior object")),
        proposal(schema.addField<int>("proposal", "archive ID for distribution used to draw samples")),
        dataSquaredNorm(schema.addField<double>("joint.r", "squared norm of weighted data vector")),
        parameterDefinition(schema.addField<int>("parameterdefinition", "archive ID for ParameterDefinition"))
    {
        schema.getCitizen().markPersistent();
    }

};

} // anonymous

// factory class used to unpersist a SampleSet
class SampleSetFactory : public tbl::io::PersistableFactory {
public:

    virtual PTR(tbl::io::Persistable)
    read(InputArchive const & archive, CatalogVector const & catalogs) const {
        LSST_ARCHIVE_ASSERT(catalogs.size() == 2u);
        LSST_ARCHIVE_ASSERT(catalogs.back().size() == 1u);
        SampleSetPersistenceKeys const & keys2 = SampleSetPersistenceKeys::get();
        LSST_ARCHIVE_ASSERT(catalogs.back().getSchema() == keys2.schema);
        tbl::BaseRecord const & record2 = catalogs.back().front();
        PTR(SampleSet) result(
            new SampleSet(
                archive.get<ParameterDefinition>(record2.get(keys2.parameterDefinition)),
                catalogs.front()
            )
        );
        result->setDataSquaredNorm(record2.get(keys2.dataSquaredNorm));
        result->setProposal(archive.get<MixtureBase>(record2.get(keys2.proposal)));
        result->_prior = archive.get<Prior>(record2.get(keys2.prior));
        return result;
    }

    explicit SampleSetFactory(std::string const & name) : tbl::io::PersistableFactory(name) {}
};

static std::string getSampleSetPersistenceName() { return "SampleSet"; }

// constructor for this instance registers the factor in a singleton in afw::table::io
static SampleSetFactory registration(getSampleSetPersistenceName());

std::string SampleSet::getPersistenceName() const {
    return getSampleSetPersistenceName();
}

std::string SampleSet::getPythonModule() const {
    return "lsst.meas.multifit";
}

void SampleSet::write(OutputArchiveHandle & handle) const {
    SampleSetPersistenceKeys const & keys2 = SampleSetPersistenceKeys::get();
    tbl::BaseCatalog catalog1 = handle.makeCatalog(_keys.schema);
    catalog1.insert(catalog1.end(), _records.begin(), _records.end());
    handle.saveCatalog(catalog1);
    tbl::BaseCatalog catalog2 = handle.makeCatalog(keys2.schema);
    PTR(tbl::BaseRecord) record2 = catalog2.addNew();
    record2->set(keys2.parameterDefinition, handle.put(_parameterDefinition));
    record2->set(keys2.dataSquaredNorm, _dataSquaredNorm);
    record2->set(keys2.proposal, handle.put(_proposal));
    record2->set(keys2.prior, handle.put(_prior));
    handle.saveCatalog(catalog2);
}

}}} // namespace lsst::meas::multifit
