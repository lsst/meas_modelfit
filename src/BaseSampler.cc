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

#include "lsst/pex/exceptions.h"
#include "lsst/afw/table/io/CatalogVector.h"
#include "lsst/afw/table/io/InputArchive.h"
#include "lsst/afw/table/io/OutputArchive.h"
#include "lsst/meas/multifit/BaseSampler.h"
#include "lsst/meas/multifit/priors.h"

namespace lsst { namespace meas { namespace multifit {

SamplePoint::SamplePoint(int nonlinearDim, int linearDim) :
    joint(linearDim),
    marginal(std::numeric_limits<Pixel>::quiet_NaN()),
    proposal(1.0),
    parameters(Vector::Zero(nonlinearDim))
{}

SampleSet::SampleSet(int nonlinearDim, int linearDim) :
    _nonlinearDim(nonlinearDim), _linearDim(linearDim)
{}


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
    _samples.push_back(p);
    if (_prior) {
        _prior->apply(_samples.back().joint, _samples.back().parameters);
    }
}

void SampleSet::applyPrior(PTR(Prior) const & prior) {
    for (iterator i = begin(); i != end(); ++i) {
        i->marginal = prior->apply(i->joint, i->parameters);
    }
    _prior = prior;
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
            "Cannot compute expectation values without marginalizing over amplitudes"
        );
    }
    // n.b. we could probably write this as a single-pass algorithm, which would only involve
    // evaluating the expectation functor once for each point, but the single-pass algorithm
    // is a lot easier to read than a robust single-pass algorithm.
    double alpha = 0.0;
    Eigen::VectorXd r = Eigen::VectorXd::Zero(functor.outputDim);
    for (const_iterator i = begin(); i != end(); ++i) {
        double w = i->marginal / i->proposal;
        alpha += w;
        r += w * functor(*i, *_prior).cast<double>();
    }
    r /= alpha;
    if (mcCov) {
        mcCov->setZero();
        alpha /= size();
        double alphaVar = 0.0;
        for (const_iterator i = begin(); i != end(); ++i) {
            double w = i->marginal / i->proposal;
            alphaVar += (w - alpha);
            Eigen::VectorXd delta = functor(*i, *_prior).cast<double>() - r;
            mcCov->selfadjointView<Eigen::Lower>().rankUpdate(delta, 1.0);
        }
        mcCov->selfadjointView<Eigen::Lower>().rankUpdate(r, alphaVar);
        *mcCov = mcCov->selfadjointView<Eigen::Lower>();
        *mcCov /= (alpha * alpha * size() * (size() - 1));
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

namespace tbl = afw::table;

class SampleSetPersistenceHelper {
public:
    tbl::Schema schema;
    tbl::Key<Pixel> joint_r;
    tbl::Key< tbl::Array<Pixel> > joint_mu;
    tbl::Key< tbl::Covariance<Pixel> > joint_fisher;
    tbl::Key<Pixel> marginal;
    tbl::Key<Pixel> proposal;
    tbl::Key< tbl::Array<Pixel> > parameters;

    SampleSetPersistenceHelper(int nonlinearDim, int linearDim) :
        schema(),
        joint_r(schema.addField<Pixel>("joint.r", "1/2 chi^2 at maximum likelihood point")),
        joint_mu(schema.addField< tbl::Array<Pixel> >("joint.mu", "maximum likelihood amplitude vector",
                                                      linearDim)),
        joint_fisher(schema.addField< tbl::Covariance<Pixel> >("joint.fisher", "amplitude Fisher matrix",
                                                               linearDim)),
        marginal(schema.addField<Pixel>("marginal", "marginal posterior value at sample point")),
        proposal(schema.addField<Pixel>("proposal", "density of the distribution used to draw samples")),
        parameters(schema.addField< tbl::Array<Pixel> >("parameters", "nonlinear parameters at this point",
                                                        nonlinearDim))
    {}

    explicit SampleSetPersistenceHelper(tbl::Schema const & schema_) :
        schema(schema_),
        joint_r(schema["joint.r"]),
        joint_mu(schema["joint.mu"]),
        joint_fisher(schema["joint.fisher"]),
        marginal(schema["marginal"]),
        proposal(schema["proposal"]),
        parameters(schema["parameters"])
    {}

};

class SampleSetFactory : public tbl::io::PersistableFactory {
public:

    virtual PTR(tbl::io::Persistable)
    read(InputArchive const & archive, CatalogVector const & catalogs) const {
        LSST_ARCHIVE_ASSERT(catalogs.size() == 1u);
        SampleSetPersistenceHelper const keys(catalogs.front().getSchema());
        int const linearDim = keys.joint_mu.getSize();
        int const nonlinearDim = keys.parameters.getSize();
        PTR(SampleSet) result(new SampleSet(nonlinearDim, linearDim));
        result->reserve(catalogs.front().size());
        for (
            tbl::BaseCatalog::const_iterator i = catalogs.front().begin();
            i != catalogs.front().end();
            ++i
        ) {
            SamplePoint p(nonlinearDim, linearDim);
            p.joint.r = i->get(keys.joint_r);
            p.joint.mu = i->get(keys.joint_mu).asEigen();
            p.joint.fisher = i->get(keys.joint_fisher);
            p.marginal = i->get(keys.marginal);
            p.proposal = i->get(keys.proposal);
            p.parameters = i->get(keys.parameters).asEigen();
            result->add(p);
        }
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
    SampleSetPersistenceHelper const keys(_nonlinearDim, _linearDim);
    tbl::BaseCatalog catalog = handle.makeCatalog(keys.schema);
    catalog.reserve(size());
    for (const_iterator i = begin(); i != end(); ++i) {
        PTR(tbl::BaseRecord) record = catalog.addNew();
        record->set(keys.joint_r, i->joint.r);
        (*record)[keys.joint_mu].asEigen() = i->joint.mu;
        record->set(keys.joint_fisher, i->joint.fisher);
        record->set(keys.marginal, i->marginal);
        record->set(keys.proposal, i->proposal);
        (*record)[keys.parameters].asEigen() = i->parameters;
    }
    handle.saveCatalog(catalog);
}

}}} // namespace lsst::meas::multifit
