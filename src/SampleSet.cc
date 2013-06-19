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
    jointFisher(schema.addField<samples::ArrayTag>("joint.fisher", "amplitude Fisher matrix", linearDim)),
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
    assert(*record.getElement(jointFisher) == joint.fisher(0,0));
}

SampleSet::SampleSet(int nonlinearDim, int linearDim, std::string const & ellipseType) :
    _keys(nonlinearDim, linearDim),
    _records(_keys.schema),
    _dataSquaredNorm(0.0),
    _ellipseType(ellipseType),
    _prior()
{}

SampleSet::SampleSet(tbl::BaseCatalog const & records, std::string const & ellipseType) :
    _keys(records.getSchema()),
    _records(records),
    _dataSquaredNorm(0.0),
    _ellipseType(ellipseType),
    _prior()
{}

void SampleSet::setEllipseType(std::string const & ellipseType) {
    if (ellipseType == _ellipseType) return;
    PTR(afw::geom::ellipses::BaseCore) oldEllipse = afw::geom::ellipses::BaseCore::make(_ellipseType);
    PTR(afw::geom::ellipses::BaseCore) newEllipse = afw::geom::ellipses::BaseCore::make(ellipseType);
    ndarray::Array<samples::Scalar,2,1> parameters = _records.getColumnView()[_keys.parameters];
    for (ndarray::Array<samples::Scalar,2,1>::Iterator i = parameters.begin(); i != parameters.end(); ++i) {
        oldEllipse->setParameterVector(i->asEigen().segment<3>(0));
        *newEllipse = *oldEllipse;
        i->asEigen().segment<3>(0) = newEllipse->getParameterVector();
    }
    _ellipseType = ellipseType;
}

afw::geom::ellipses::Ellipse SampleSet::interpret(
    samples::Vector const & parameters,
    afw::geom::Point2D const & center
) const {
    afw::geom::ellipses::Ellipse result(
        afw::geom::ellipses::BaseCore::make(_ellipseType, parameters.segment<3>(0)),
        center
    );
    if (parameters.size() >= 5) {
        result.setCenter(afw::geom::Point2D(parameters.segment<2>(3)));
    }
    return result;
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

double SampleSet::applyPrior(PTR(Prior) const & prior) {
    for (tbl::BaseCatalog::iterator i = _records.begin(); i != _records.end(); ++i) {
        i->set(_keys.marginal, prior->apply(_keys.getJoint(*i), _keys.getParameters(*i)));
        // for numerical reasons, in the first pass, we set w_i = ln(m_i/q_i);
        // note that i->proposal == -ln(q_i) and i->marginal == -ln(m_i)
        i->set(_keys.weight, i->get(_keys.proposal) - i->get(_keys.marginal));
    }
    // sort by ascending probability, so when we accumulate, we add small numbers together
    // before adding them to large numbers
    _records.sort(_keys.weight);
    // we now compute z, the arithmetic mean of ln(m_i/q_i)
    double z = 0.0;
    for (tbl::BaseCatalog::iterator i = _records.begin(); i != _records.end(); ++i) {
        z += i->get(_keys.weight);
    }
    z /= _records.size();
    // we now subtract z from w_i and exponentiate, accumulating the sums
    // this now makes w_i = e^{-z} m_i/q_i, which is proportional to the
    // desired m_i/q_i.
    double wSum = 0.0;
    for (tbl::BaseCatalog::iterator i = _records.begin(); i != _records.end(); ++i) {
        wSum += (*i)[_keys.weight] = std::exp(i->get(_keys.weight) - z);
    }
    // finally, we normalize w_i...
    for (tbl::BaseCatalog::iterator i = _records.begin(); i != _records.end(); ++i) {
        (*i)[_keys.weight] /= wSum;
    }
    _prior = prior;
    // ..and return the log of wSum, corrected for the z term we took out earlier,
    // and including the r/2 term we've ignored all along.
    return 0.5*_dataSquaredNorm - z - std::log(wSum / _records.size());
}

void SampleSet::dropPrior() {
    for (tbl::BaseCatalog::iterator i = _records.begin(); i != _records.end(); ++i) {
        (*i)[_keys.weight] = (*i)[_keys.marginal] = std::numeric_limits<double>::quiet_NaN();
    }
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

// schema and keys for single-record (per SampleSet), dimension-independent catalog
class SampleSetPersistenceKeys : private boost::noncopyable {
public:
    tbl::Schema schema;
    tbl::Key<double> dataSquaredNorm;
    tbl::Key<std::string> ellipseType;

    static SampleSetPersistenceKeys const & get() {
        static SampleSetPersistenceKeys const instance;
        return instance;
    }

private:

    SampleSetPersistenceKeys() :
        schema(),
        dataSquaredNorm(schema.addField<double>("joint.r", "squared norm of weighted data vector")),
        ellipseType(schema.addField<std::string>("ellipsetype", "name of ellipse parametrization", 48))
    {
        schema.getCitizen().markPersistent();
    }

};

class SampleSetFactory : public tbl::io::PersistableFactory {
public:

    virtual PTR(tbl::io::Persistable)
    read(InputArchive const & archive, CatalogVector const & catalogs) const {
        LSST_ARCHIVE_ASSERT(catalogs.size() == 2u);
        LSST_ARCHIVE_ASSERT(catalogs.back().size() == 1u);
        SampleSetPersistenceKeys const & keys2 = SampleSetPersistenceKeys::get();
        LSST_ARCHIVE_ASSERT(catalogs.back().getSchema() == keys2.schema);
        tbl::BaseRecord const & record2 = catalogs.back().front();
        PTR(SampleSet) result(new SampleSet(catalogs.front(), record2.get(keys2.ellipseType)));
        result->setDataSquaredNorm(record2.get(keys2.dataSquaredNorm));
        return result;
    }

    explicit SampleSetFactory(std::string const & name) : tbl::io::PersistableFactory(name) {}
};

std::string getSampleSetPersistenceName() { return "SampleSet"; }

SampleSetFactory registration(getSampleSetPersistenceName());

} // anonymous

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
    record2->set(keys2.ellipseType, _ellipseType);
    record2->set(keys2.dataSquaredNorm, _dataSquaredNorm);
    handle.saveCatalog(catalog2);
}

}}} // namespace lsst::meas::multifit
