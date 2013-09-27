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

#ifndef LSST_MEAS_MULTIFIT_SampleSet_h_INCLUDED
#define LSST_MEAS_MULTIFIT_SampleSet_h_INCLUDED

#include "lsst/afw/table/BaseRecord.h"
#include "lsst/afw/table/BaseTable.h"
#include "lsst/afw/table/Catalog.h"
#include "lsst/afw/table/io/Persistable.h"
#include "lsst/meas/multifit/constants.h"
#include "lsst/meas/multifit/LogGaussian.h"
#include "lsst/meas/multifit/priors.h"
#include "lsst/meas/multifit/parameters.h"
#include "lsst/meas/multifit/Mixture.h"
#include "lsst/meas/multifit/KernelDensityEstimator.h"

namespace lsst { namespace meas { namespace multifit {

class ExpectationFunctor;

class SampleSetKeys {
public:
    afw::table::Schema schema;
    samples::ArrayKey jointGrad;
    samples::ArrayKey jointFisher;
    samples::ScalarKey marginal;
    samples::ScalarKey proposal;
    samples::ScalarKey weight;
    samples::ArrayKey parameters;

    /// Return the number of nonlinear parameters
    int getParameterDim() const { return parameters.getSize(); }

    /// Return the number of linear coefficients
    int getCoefficientDim() const { return jointGrad.getSize(); }

    /// Extract a LogGaussian object from the given record
    LogGaussian getJoint(afw::table::BaseRecord const & record) const;

    /// Set LogGaussian fields in the given record
    void setJoint(afw::table::BaseRecord & record, LogGaussian const & joint);

#ifndef SWIG

    /// Extract the nonlinear parameter vector as a const Eigen Map object
    samples::VectorCMap getParameters(afw::table::BaseRecord const & record) const {
        return samples::VectorCMap(record.getElement(parameters), getParameterDim());
    }

    /// Set nonlinear parameter fields in the given record
    void setParameters(afw::table::BaseRecord & record, samples::Vector const & parameters_) const {
        samples::VectorMap(record.getElement(parameters), getParameterDim()) = parameters_;
    }

#endif // !SWIG

    SampleSetKeys(int parameterDim, int coefficientDim);

    explicit SampleSetKeys(afw::table::Schema const & schema_);
};

/**
 *  @brief Representation of a probability distribution as a set of Monte Carlo samples that
 *         distinguishes linear coefficients from other nonlinear parameters.
 *
 *  For linear coefficients @f$\alpha@f$ and nonlinear parameters @f$\theta@f$,
 *  each element @f$n@f$ in a SampleSet contains:
 *   - the nonlinear parameters @f$\theta_n@f$ at that point
 *   - the joint likelihood @f$P(D|\alpha,\theta_n) = e^{-L_n(\alpha)}@f$ (see LogGaussian)
 *   - the nonnormalized marginal posterior
 *     @f$m_n=P(\theta_n|D)P(D)=P(D|\theta_n)P(\theta_n)@f$ (the Bayesian evidence @f$P(D)@f$
 *     is the normalization factor), as obtained applying a Prior to the joint
 *     likelihood at each SamplePoint.
 *   - the density @f$q_n@f$ of the distribution from which the samples were drawn
 *
 *  Together with the Prior, these indirectly provide a full representation of the
 *  joint posterior (see ExpectationFunctor for mathematical details of how expectation values
 *  on the full posterior are computed).
 */
class SampleSet : public afw::table::io::PersistableFacade<SampleSet>, public afw::table::io::Persistable {
public:

    /**
     *  @brief Initialize the SampleSet with the given parameter dimensions.
     *
     *  Any SamplePoints added to the SampleSet must have the same dimensions.
     *
     *  @param[in] parameterDefinition    Object that defines the interpretation of the nonlinear parameters.
     *  @param[in] coefficientDim      Number of coefficient coefficients (dimension of SamplePoint::joint).
     */
    SampleSet(PTR(ParameterDefinition const) parameterDefinition, int coefficientDim);

    /**
     *  @brief Construct a SampleSet from an existing catalog of sample records.
     *
     *  @param[in] parameterDefinition    Object that defines the interpretation of the nonlinear parameters.
     *  @param[in] records        Catalog of records with schema compatible with SampleSetKeys
     */
    SampleSet(PTR(ParameterDefinition const) parameterDefinition, afw::table::BaseCatalog const & records);

    /// Return the number of nonlinear parameters
    int getParameterDim() const { return _keys.getParameterDim(); }

    /// Return the number of linear coefficients
    int getCoefficientDim() const { return _keys.getCoefficientDim(); }

    /// Return the object that defines how to interpret the nonlinear parameters
    PTR(ParameterDefinition const) getParameterDefinition() const { return _parameterDefinition; }

    /// Convert the nonlinear parameters according to a new ParameterDefinition
    void setParameterDefinition(PTR(ParameterDefinition const) parameterDefinition);

    /**
     *  @brief Return the squared norm of weighted data values
     *
     *  This is the @f$r@f$ quantity logically associated with each LogGaussian, but
     *  because it is a constant shared by all samples, it is not stored with each
     *  LogGaussian.
     */
    double getDataSquaredNorm() const { return _dataSquaredNorm; }

    /**
     *  @brief Set the squared norm of weighted data values
     *
     *  @sa getDataSquaredNorm()
     */
    void setDataSquaredNorm(double r) { _dataSquaredNorm = r; }

    /**
     *  @brief Return the distribution from which the nonlinear parameters were originally drawn.
     */
    PTR(MixtureBase const) getProposal() const { return _proposal; }

    /**
     *  @brief Set the distribution from which the nonlinear parameters were originally drawn.
     *
     *  This does not recompute the 'proposal' values for each sample.
     */
    void setProposal(PTR(MixtureBase const) proposal) { _proposal = proposal; }

    /// @brief Return an afw::table::BaseCatalog representation of the SampleSet.
    afw::table::BaseCatalog getCatalog() const { return _records; }

    /**
     *  @brief Create a 1-d estimate of the density of a single parameter (marginalized over the others)
     */
    ndarray::Array<double,1,1> computeDensity(KernelDensityEstimatorControl const & ctrl) const;

    /**
     *  @brief Create a 2-d estimate of the density of two parameters (marginalized over the others)
     *
     *  Note that the returned object follows the usual [Y][X] convention for ndarray images;
     *  the range specified in ctrlY defines the rows (zeroth dimension) of the returned Array, while
     *  ctrlX specifies the columns (first dimension).
     */
    ndarray::Array<double,2,2> computeDensity(
        KernelDensityEstimatorControl const & ctrlX,
        KernelDensityEstimatorControl const & ctrlY
    ) const;

    /**
     *  @brief Compute the normalized perplexity of the samples
     *
     *  This is the exponential of the Shannon entropy of the weights, divided by the number of samples:
     *  @f[
     *     \frac{1}{N}\exp\left(-\sum_i^N w_i \ln w_i\right)
     *  @f]
     */
    double computeNormalizedPerplexity() const;

    /**
     *  @brief Compute the fraction of the effective sample size over the actual sample size
     *
     *  The effective sample size is the inverse of the sum of squares of weights:
     *  @f[
     *     \frac{1}{N}\left(\sum_i^N w_i^2 \right)^{-1}
     *  @f]
     */
    double computeEffectiveSampleSizeFraction() const;

    /// Return the number of samples.
    std::size_t size() const { return _records.size(); }

    /**
     *  @brief Add a new sample point to the SampleSet.
     *
     *  If a prior has already been applied to the SampleSet, new points may not
     *  be added (throws LogicErrorException).
     */
    void add(LogGaussian const & joint, samples::Scalar proposal, samples::Vector const & parameters);

    /**
     *  @brief Attach the given prior to the SampleSet and apply it to all existing samples.
     *
     *  Attaching a prior recomputes the "marginal" and "weight" values for each SamplePoint,
     *  sorts the samples by weight, and prevents new SamplePoints from being added.
     *
     *  Return the negative log of the sum of weights before normalizing them, which is the
     *  negative log of the Bayesian "evidence" for a normalized prior.
     *
     *  The clip parameter is used to clip samples with weights less than clip*max(weights);
     *  if set to a reasonable value, this can dramatically reduce the number of samples
     *  that need to be saved without affecting the distribution when using a naive sampler.
     *  If the given value is less than std::numeric_limits<double>::min(), it will be reset
     *  to that value to avoid numerical issues.
     */
    double applyPrior(PTR(Prior const) prior, double clip=0);

    /// Remove the prior from the SampleSet, allowing new SamplePoints to be added.
    void dropPrior();

    /// Remove all records and drop the prior.
    void clear();

    /**
     *  @brief Compute quantiles of a single nonlinear parameters
     *
     *  A quantile is the point at which the cumulative distribution reaches a certain value.
     *  For instance, the median is the quantile at 0.5, and a pair of quantiles at (0.05, 0.95)
     *  is a 90% confidence interval.
     *
     *  @param[in] fractions       A sorted array of fractions (floating point numbers in the range [0,1])
     *                             at which to compute a quantile value (e.g. 0.5 == median).
     *  @param[in] parameterIndex  Index of the parameter over which to compute the quantile; an
     *                             integer between 0 and (getParameterDim()-1).  All other parameters
     *                             are ignored (meaning they areeffectively marginalized over).
     */
    samples::Vector computeQuantiles(samples::Vector const & fractions, int parameterIndex) const;

    /**
     *  @brief Compute the same quantiles for all nonlinear parameters
     *
     *  A quantile is the point at which the cumulative distribution reaches a certain value.
     *  For instance, the median is the quantile at 0.5, and a pair of quantiles at (0.05, 0.95)
     *  is a 90% confidence interval.
     *
     *  @param[in] fractions       A sorted array of fractions (floating point numbers in the range [0,1])
     *                             at which to compute a quantile value (e.g. 0.5 == median).
     */
    samples::Matrix computeQuantiles(samples::Vector const & fractions) const;

    /**
     *  @brief Compute an expectation integral and optionally its Monte Carlo covariance.
     *
     *  See ExpectationFunctor for details of the calculation.
     *
     *  If the Monte Carlo covariance matrix is requested, it will represent the uncertainty due only to
     *  the finite number of samples and non-optimality of the proposal distribution, not the uncertainty
     *  due to the width of the distribution itself.
     *
     *  A prior must be attached before computeExpecation is called.
     *
     *  See ExpectationFunctor for more information.
     */
    samples::Vector computeExpectation(ExpectationFunctor const & functor, samples::Matrix * mcCov=0) const;

    /**
     *  @brief Compute the empirical mean of the marginal distribution.
     *
     *  Using the notation of ExpectationFunctor, this is the expectation value for
     *  @f$f(\alpha,\theta)=\theta@f$.
     *
     *  As with computeExpectation, the optional "mcCov" output represents only the uncertainty due to the
     *  finite number of samples, and is not the covariance matrix of the distribution.
     *
     *  A prior must be attached before computeExpecation is called.
     */
    samples::Vector computeMean(samples::Matrix * mcCov=0) const;

    /**
     *  @brief Compute the empirical covariance of the marginal distribution.
     *
     *  Using the notation of ExpectationFunctor, this is the expectation value for
     *  @f$f(\alpha,\theta)=\theta\theta^T@f$.
     *
     *  A prior must be attached before computeExpecation is called.
     */
    samples::Matrix computeCovariance(samples::Vector const & mean) const;

    /// @copydoc computeCovariance
    samples::Matrix computeCovariance() const { return computeCovariance(computeMean()); }

    bool isPersistable() const { return true; }

protected:

    friend class SampleSetFactory;

    virtual std::string getPersistenceName() const;

    virtual std::string getPythonModule() const;

    virtual void write(OutputArchiveHandle & handle) const;

private:
    SampleSetKeys _keys;
    afw::table::BaseCatalog _records;
    double _dataSquaredNorm;
    PTR(ParameterDefinition const) _parameterDefinition;
    PTR(Prior const) _prior;
    PTR(MixtureBase const) _proposal;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_SampleSet_h_INCLUDED
