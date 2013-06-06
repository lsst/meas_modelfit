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
#include "lsst/afw/table/io/CatalogVector.h"
#include "lsst/afw/table/io/InputArchive.h"
#include "lsst/afw/table/io/OutputArchive.h"
#include "lsst/meas/multifit/SampleSet.h"
#include "lsst/meas/multifit/ExpectationFunctor.h"
#include "lsst/meas/multifit/priors.h"

namespace tbl = lsst::afw::table;

namespace lsst { namespace meas { namespace multifit {

SamplePoint::SamplePoint(int nonlinearDim, int linearDim) :
    joint(linearDim),
    marginal(std::numeric_limits<Pixel>::quiet_NaN()),
    proposal(0.0),
    weight(std::numeric_limits<Pixel>::quiet_NaN()),
    parameters(Vector::Zero(nonlinearDim))
{}

SampleSet::SampleSet(int nonlinearDim, int linearDim, std::string const & ellipseType) :
    _nonlinearDim(nonlinearDim), _linearDim(linearDim), _ellipseType(ellipseType)
{}

void SampleSet::setEllipseType(std::string const & ellipseType) {
    if (ellipseType == _ellipseType) return;
    PTR(afw::geom::ellipses::BaseCore) oldEllipse = afw::geom::ellipses::BaseCore::make(_ellipseType);
    PTR(afw::geom::ellipses::BaseCore) newEllipse = afw::geom::ellipses::BaseCore::make(ellipseType);
    for (iterator i = begin(); i != end(); ++i) {
        oldEllipse->setParameterVector(i->parameters.segment<3>(0).cast<double>());
        *newEllipse = *oldEllipse;
        i->parameters.segment<3>(0) = newEllipse->getParameterVector().cast<float>();
    }
    _ellipseType = ellipseType;
}

afw::geom::ellipses::Ellipse SampleSet::interpret(
    Eigen::VectorXd const & parameters,
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
    ndarray::Array<double,1,1> result = ndarray::allocate(ctrl.getCount());
    result.deep() = 0.0;
    for (const_iterator s = begin(); s != end(); ++s) {
        ctrl.apply(result, s->parameters[ctrl.getIndex()], s->weight);
    }
    return result;
}

ndarray::Array<double,2,2> SampleSet::computeDensity(
    KernelDensityEstimatorControl const & ctrlX,
    KernelDensityEstimatorControl const & ctrlY
) const {
    ndarray::Array<double,2,2> result = ndarray::allocate(ctrlY.getCount(), ctrlX.getCount());
    ndarray::Array<double,1,1> tmpY = ndarray::allocate(ctrlY.getCount());
    ndarray::Array<double,1,1> tmpX = ndarray::allocate(ctrlX.getCount());
    result.deep() = 0.0;
    for (const_iterator s = begin(); s != end(); ++s) {
        tmpY.deep() = 0.0;
        tmpX.deep() = 0.0;
        ctrlY.apply(tmpY, s->parameters[ctrlY.getIndex()], 1.0);
        ctrlX.apply(tmpX, s->parameters[ctrlX.getIndex()], 1.0);
        result.asEigen() += tmpY.asEigen() * tmpX.asEigen().adjoint();
    }
    return result;
}

void SampleSet::add(SamplePoint const & p) {
    if (p.parameters.size() != _nonlinearDim) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            (boost::format("Incorrect nonlinear dimension for SamplePoint: expected %d, got %d")
             % _nonlinearDim % p.parameters.size()).str()
        );
    }
    if (p.joint.mu.size() != _linearDim) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            (boost::format("Incorrect linear dimension for SamplePoint: expected %d, got %d")
             % _nonlinearDim % p.joint.mu.size()).str()
        );
    }
    if (_prior) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "Cannot add new samples after applyPrior() has been called; use dropPrior() to revert"
        );
    }
    _samples.push_back(p);
}

double SampleSet::applyPrior(PTR(Prior) const & prior) {
    for (iterator i = begin(); i != end(); ++i) {
        i->marginal = prior->apply(i->joint, i->parameters);
        // for numerical reasons, in the first pass, we set w_i = ln(m_i/q_i);
        // note that i->proposal == -ln(q_i) and i->marginal == -ln(m_i)
        i->weight = i->proposal - i->marginal;
    }
    // sort by ascending probability, so when we accumulate, we add small numbers together
    // before adding them to large numbers
    _samples.sort();
    // we now compute z, the arithmetic mean of ln(m_i/q_i)
    double z = 0.0;
    for (iterator i = begin(); i != end(); ++i) {
        z += i->weight;
    }
    z /= _samples.size();
    // we now subtract z from w_i and exponentiate, accumulating the sums
    // this now makes w_i = e^{-z} m_i/q_i, which is proportional to the
    // desired m_i/q_i.
    double wSum = 0.0;
    for (iterator i = begin(); i != end(); ++i) {
        wSum += i->weight = std::exp(i->weight - z);
    }
    // finally, we normalize w_i...
    for (iterator i = begin(); i != end(); ++i) {
        i->weight /= wSum;
    }
    _prior = prior;
    // ..and return the log of wSum, corrected for the z term we took out earlier
    return -z - std::log(wSum / _samples.size());
}

void SampleSet::dropPrior() {
    for (iterator i = begin(); i != end(); ++i) {
        i->weight = i->marginal = std::numeric_limits<double>::quiet_NaN();
    }
    _prior.reset();
}

namespace {

class MeanExpectationFunctor : public ExpectationFunctor {
public:
    explicit MeanExpectationFunctor(int dim_) : ExpectationFunctor(dim_) {}

    virtual Eigen::VectorXd operator()(SamplePoint const & sample, Prior const & prior) const {
        return sample.parameters.cast<double>();
    }
};

class CovarianceExpectationFunctor : public ExpectationFunctor {
public:
    explicit CovarianceExpectationFunctor(Eigen::VectorXd const & mean) :
        ExpectationFunctor(mean.size() * (mean.size() - 1)), _mean(mean)
    {}

    virtual Eigen::VectorXd operator()(SamplePoint const & sample, Prior const & prior) const {
        Eigen::VectorXd r = Eigen::VectorXd::Zero(outputDim);
        int n = 0;
        for (int i = 0; i < _mean.size(); ++i) {
            for (int j = 0; j <= i; ++j) {
                r[n] = (sample.parameters[i] - _mean[i]) * (sample.parameters[j] - _mean[j]);
                ++n;
            }
        }
        return r;
    }

    Eigen::MatrixXd expand(Eigen::VectorXd const & x) const {
        Eigen::MatrixXd r = Eigen::MatrixXd::Zero(_mean.size(), _mean.size());
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
    Eigen::VectorXd const & _mean;
};

} // anonymous

Eigen::VectorXd SampleSet::computeExpectation(
    ExpectationFunctor const & functor,
    Eigen::MatrixXd * mcCov
) const {
    if (!_prior) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "Cannot compute expectation values without attaching a prior"
        );
    }
    Eigen::VectorXd r = Eigen::VectorXd::Zero(functor.outputDim);
    for (const_iterator i = begin(); i != end(); ++i) {
        r += i->weight * functor(*i, *_prior);
    }
    if (mcCov) {
        mcCov->setZero();
        double w2 = 0.0;
        for (const_iterator i = begin(); i != end(); ++i) {
            Eigen::VectorXd delta = functor(*i, *_prior) - r;
            mcCov->selfadjointView<Eigen::Lower>().rankUpdate(delta, i->weight);
            w2 += i->weight * i->weight;
        }
        *mcCov = mcCov->selfadjointView<Eigen::Lower>();
        *mcCov /= (1.0 - w2);
    }
    return r;
}

Eigen::VectorXd SampleSet::computeMean(Eigen::MatrixXd * mcCov) const {
    MeanExpectationFunctor f(_nonlinearDim);
    return computeExpectation(f, mcCov);
}

Eigen::MatrixXd SampleSet::computeCovariance(Eigen::VectorXd const & mean) const {
    CovarianceExpectationFunctor f(mean);
    return f.expand(computeExpectation(f));
}

namespace {

// Schema and keys for multi-record, dimension-dependent catalog
class SampleSetPersistenceHelper1 {
public:
    tbl::Schema schema;
    tbl::Key<Pixel> joint_r;
    tbl::Key< tbl::Array<Pixel> > joint_mu;
    tbl::Key< tbl::Covariance<Pixel> > joint_fisher;
    tbl::Key<Pixel> marginal;
    tbl::Key<Pixel> proposal;
    tbl::Key<double> weight;
    tbl::Key< tbl::Array<Pixel> > parameters;

    SampleSetPersistenceHelper1(int nonlinearDim, int linearDim) :
        schema(),
        joint_r(schema.addField<Pixel>("joint.r", "1/2 chi^2 at maximum likelihood point")),
        joint_mu(schema.addField< tbl::Array<Pixel> >("joint.mu", "maximum likelihood amplitude vector",
                                                      linearDim)),
        joint_fisher(schema.addField< tbl::Covariance<Pixel> >("joint.fisher", "amplitude Fisher matrix",
                                                               linearDim)),
        marginal(schema.addField<Pixel>("marginal", "negative log marginal posterior value at sample point")),
        proposal(schema.addField<Pixel>("proposal",
                                        "negative log density of the distribution used to draw samples")),
        weight(schema.addField<double>("weight", "normalized Monte Carlo weight")),
        parameters(schema.addField< tbl::Array<Pixel> >("parameters", "nonlinear parameters at this point",
                                                        nonlinearDim))
    {}

    explicit SampleSetPersistenceHelper1(tbl::Schema const & schema_) :
        schema(schema_),
        joint_r(schema["joint.r"]),
        joint_mu(schema["joint.mu"]),
        joint_fisher(schema["joint.fisher"]),
        marginal(schema["marginal"]),
        proposal(schema["proposal"]),
        weight(schema["weight"]),
        parameters(schema["parameters"])
    {}

};

// schema and keys for single-record (per SampleSet), dimension-independent catalog
class SampleSetPersistenceHelper2 : private boost::noncopyable {
public:
    tbl::Schema schema;
    tbl::Key<std::string> ellipseType;

    static SampleSetPersistenceHelper2 const & get() {
        static SampleSetPersistenceHelper2 const instance;
        return instance;
    }

private:

    SampleSetPersistenceHelper2() :
        schema(),
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
        SampleSetPersistenceHelper1 const keys1(catalogs.front().getSchema());
        SampleSetPersistenceHelper2 const & keys2 = SampleSetPersistenceHelper2::get();
        LSST_ARCHIVE_ASSERT(catalogs.back().getSchema() == keys2.schema);
        int const linearDim = keys1.joint_mu.getSize();
        int const nonlinearDim = keys1.parameters.getSize();
        tbl::BaseRecord const & record2 = catalogs.back().front();
        PTR(SampleSet) result(new SampleSet(nonlinearDim, linearDim, record2.get(keys2.ellipseType)));
        for (
            tbl::BaseCatalog::const_iterator i = catalogs.front().begin();
            i != catalogs.front().end();
            ++i
        ) {
            SamplePoint p(nonlinearDim, linearDim);
            p.joint.r = i->get(keys1.joint_r);
            p.joint.mu = i->get(keys1.joint_mu).asEigen();
            p.joint.fisher = i->get(keys1.joint_fisher);
            p.marginal = i->get(keys1.marginal);
            p.proposal = i->get(keys1.proposal);
            p.weight = i->get(keys1.weight);
            p.parameters = i->get(keys1.parameters).asEigen();
            result->add(p);
        }
        return result;
    }

    explicit SampleSetFactory(std::string const & name) : tbl::io::PersistableFactory(name) {}
};

std::string getSampleSetPersistenceName() { return "SampleSet"; }

SampleSetFactory registration(getSampleSetPersistenceName());

void fillCatalog(
    SampleSet const & samples,
    tbl::BaseCatalog & catalog,
    SampleSetPersistenceHelper1 const & keys
) {
    catalog.reserve(samples.size());
    for (SampleSet::const_iterator i = samples.begin(); i != samples.end(); ++i) {
        PTR(tbl::BaseRecord) record = catalog.addNew();
        record->set(keys.joint_r, i->joint.r);
        (*record)[keys.joint_mu].asEigen() = i->joint.mu;
        record->set(keys.joint_fisher, i->joint.fisher);
        record->set(keys.marginal, i->marginal);
        record->set(keys.proposal, i->proposal);
        record->set(keys.weight, i->weight);
        (*record)[keys.parameters].asEigen() = i->parameters;
    }
}

} // anonymous

tbl::BaseCatalog SampleSet::asCatalog() const {
    SampleSetPersistenceHelper1 const keys(_nonlinearDim, _linearDim);
    tbl::BaseCatalog catalog(keys.schema);
    fillCatalog(*this, catalog, keys);
    return catalog;
}

std::string SampleSet::getPersistenceName() const {
    return getSampleSetPersistenceName();
}

std::string SampleSet::getPythonModule() const {
    return "lsst.meas.multifit";
}

void SampleSet::write(OutputArchiveHandle & handle) const {
    SampleSetPersistenceHelper1 const keys1(_nonlinearDim, _linearDim);
    SampleSetPersistenceHelper2 const & keys2 = SampleSetPersistenceHelper2::get();
    tbl::BaseCatalog catalog1 = handle.makeCatalog(keys1.schema);
    fillCatalog(*this, catalog1, keys1);
    handle.saveCatalog(catalog1);
    tbl::BaseCatalog catalog2 = handle.makeCatalog(keys2.schema);
    PTR(tbl::BaseRecord) record2 = catalog2.addNew();
    record2->set(keys2.ellipseType, _ellipseType);
    handle.saveCatalog(catalog2);
}

}}} // namespace lsst::meas::multifit
