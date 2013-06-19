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
#include "lsst/meas/multifit/KernelDensityEstimator.h"

namespace lsst { namespace meas { namespace multifit {

class ExpectationFunctor;

class SampleSetKeys {
public:
    afw::table::Schema schema;
    samples::ScalarKey jointR;
    samples::ArrayKey jointMu;
    samples::ArrayKey jointFisher;
    samples::ScalarKey marginal;
    samples::ScalarKey proposal;
    samples::ScalarKey weight;
    samples::ArrayKey parameters;

    /// Return the number of nonlinear parameters
    int getNonlinearDim() const { return parameters.getSize(); }

    /// Return the number of linear amplitude parameters
    int getLinearDim() const { return jointMu.getSize(); }

    /// Extract a LogGaussian object from the given record
    LogGaussian getJoint(afw::table::BaseRecord const & record) const;

    /// Set LogGaussian fields in the given record
    void setJoint(afw::table::BaseRecord & record, LogGaussian const & joint);

#ifndef SWIG

    /// Extract the nonlinear parameter vector as a const Eigen Map object
    samples::VectorCMap getParameters(afw::table::BaseRecord const & record) const {
        return samples::VectorCMap(record.getElement(parameters), getNonlinearDim());
    }

    /// Set nonlinear parameter fields in the given record
    void setParameters(afw::table::BaseRecord & record, samples::Vector const & parameters_) const {
        samples::VectorMap(record.getElement(parameters), getNonlinearDim()) = parameters_;
    }

#endif // !SWIG

    SampleSetKeys(int nonlinearDim, int linearDim);

    explicit SampleSetKeys(afw::table::Schema const & schema_);
};

/**
 *  @brief Representation of a probability distribution as a set of Monte Carlo samples that
 *         distinguishes linear amplitude parameters from other nonlinear parameters.
 *
 *  For linear amplitudes @f$\alpha@f$ and nonlinear parameters @f$\theta@f$,
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
     *  @param[in] nonlinearDim   Number of nonlinear parameters (size of SamplePoint::parameters).
     *  @param[in] linearDim      Number of linear amplitude parameters (dimension of SamplePoint::joint).
     *  @param[in] ellipseType    Name of the afw::geom::ellipses::BaseCore subclass that defines the
     *                            ellipse part of the nonlinear parameters.
     */
    SampleSet(int nonlinearDim, int linearDim, std::string const & ellipseType);

    /**
     *  @brief Construct a SampleSet from an existing catalog of sample records.
     *
     *  @param[in] records        Catalog of records with schema compatible with SampleSetKeys
     *  @param[in] ellipseType    Name of the afw::geom::ellipses::BaseCore subclass that defines the
     *                            ellipse part of the nonlinear parameters.
     */
    SampleSet(afw::table::BaseCatalog const & records, std::string const & ellipseType);

    /// Return the number of nonlinear parameters
    int getNonlinearDim() const { return _keys.getNonlinearDim(); }

    /// Return the number of linear amplitude parameters
    int getLinearDim() const { return _keys.getLinearDim(); }

    /// Return the name of the afw::geom::ellipses type used to interpret the nonlinear parameters.
    std::string const & getEllipseType() const { return _ellipseType; }

    /**
     *  @brief Convert the ellipse part of the nonlinear parameters to a parametrization
     *
     *  @param[in] ellipseType    Name of the afw::geom::ellipses::BaseCore subclass will define the
     *                            new parametrization for the ellipse part of the nonlinear parameters.
     */
    void setEllipseType(std::string const & ellipseType);

    /**
     *  @brief Given a nonlinear parameter vector, return an Ellipse object
     *
     *  The first three parameters are always intrepreted as an afw::geom::ellipses::BaseCore object
     *  with type set by getEllipseType().  The next two parameters, if present, are used as the ellipse
     *  center; if there are only 3 parameters, the center argument will be used instead.
     */
    afw::geom::ellipses::Ellipse interpret(
        samples::Vector const & parameters,
        afw::geom::Point2D const & center=afw::geom::Point2D()
    ) const;

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
     */
    double applyPrior(PTR(Prior) const & prior);

    /// Remove the prior from the SampleSet, allowing new SamplePoints to be added.
    void dropPrior();

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
     *                             integer between 0 and (getNonlinearDim()-1).  All other parameters
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
    std::string _ellipseType;
    PTR(Prior) _prior;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_SampleSet_h_INCLUDED
