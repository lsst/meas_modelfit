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

#include <vector>

#include "lsst/afw/table/io/Persistable.h"
#include "lsst/meas/multifit/constants.h"
#include "lsst/meas/multifit/LogGaussian.h"
#include "lsst/meas/multifit/priors.h"

namespace lsst { namespace meas { namespace multifit {

class ExpectationFunctor;

/*
 *  @brief Point in a Monte-Carlo SampleSet.
 *
 *  See SampleSet for more information.
 */
class SamplePoint {
public:

    LogGaussian joint; ///< Log likelihood w.r.t. linear amplitudes @f$L_n(\alpha)@f$
    Pixel marginal;    ///< Marginal nonnormalized posterior @f$m_n@f$
    Pixel proposal;    ///< Density @f$q_n@f$ of the distribution from which the samples were drawn
    Vector parameters; ///< Nonlinear parameters @f$\theta_n@f$ at this point

    /// Initialize to zeros with the given dimensions
    SamplePoint(int nonlinearDim, int linearDim);

};

/**
 *  @brief Representation of a probability distribution as a set of Monte Carlo samples that
 *         distinguishes linear amplitude parameters from other nonlinear parameters.
 *
 *  For linear amplitudes @f$\alpha@f$ and nonlinear parameters @f$\theta@f$,
 *  each SamplePoint @f$n@f$ in a SampleSet contains:
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
    typedef std::vector<SamplePoint> Container;
public:

    typedef Container::iterator iterator;
    typedef Container::const_iterator const_iterator;

    /**
     *  @brief Initialize the SampleSet with the given parameter dimensions.
     *
     *  Any SamplePoints added to the SampleSet must have the same dimensions.
     */
    SampleSet(int nonlinearDim, int linearDim);

    //@{
    /**
     *  Iterate over SamplePoints.
     *
     *  Iterators are std::vector-based, and may be invalidated when adding new points
     *  unless the new size is less than capacity().
     */
    iterator begin() { return _samples.begin(); }
    iterator end() { return _samples.end(); }
    const_iterator begin() const { return _samples.begin(); }
    const_iterator end() const { return _samples.end(); }
    //@}

    /// Return the number of samples.
    std::size_t size() const { return _samples.size(); }

    /// Reserve space for the given total number of samples.
    void reserve(std::size_t capacity) { _samples.reserve(capacity); }

    /// Return the total number of samples that space has been allocated for.
    std::size_t capacity() const { return _samples.capacity(); }

    /**
     *  @brief Add a new sample point to the SampleSet, applying the prior
     *         if applyPrior() has been called.
     */
    void add(SamplePoint const & p);

    /**
     *  @brief Attach the given prior to the SampleSet and apply it to all existing samples.
     *
     *  Attaching a prior recomputes the "marginal" value for each SamplePoint, and causes any
     *  the prior to be applied to any future samples automatically.
     */
    void applyPrior(PTR(Prior) const & prior);

    /**
     *  @brief Compute an expectation integral and optionally its Monte Carlo covariance.
     *
     *  See ExpectationFunctor for details of the calculation.
     *
     *  If the Monte Carlo covariance matrix is requested, it will represent the uncertainty due only to
     *  the finite number of samples and non-optimality of the proposal distribution, not the uncertainty
     *  due to the width of the distribution itself.
     *
     *  See ExpectationFunctor for more information.
     */
    Eigen::VectorXd computeExpectation(ExpectationFunctor const & functor, Eigen::MatrixXd * mcCov=0) const;

    /**
     *  @brief Compute the empirical mean of the marginal distribution.
     *
     *  Using the notation of ExpectationFunctor, this is the expectation value for
     *  @f$f(\alpha,\theta)=\theta@f$.
     *
     *  As with computeExpectation, the optional "mcCov" output represents only the uncertainty due to the
     *  finite number of samples, and is not the covariance matrix of the distribution.
     */
    Eigen::VectorXd computeMean(Eigen::MatrixXd * mcCov=0) const;

    /**
     *  @brief Compute the empirical covariance of the marginal distribution.
     *
     *  Using the notation of ExpectationFunctor, this is the expectation value for
     *  @f$f(\alpha,\theta)=\theta\theta^T@f$.
     */
    Eigen::MatrixXd computeCovariance(Eigen::VectorXd const & mean) const;

    /// @copydoc computeCovariance
    Eigen::MatrixXd computeCovariance() const { return computeCovariance(computeMean()); }

    bool isPersistable() const { return true; }

protected:

    virtual std::string getPersistenceName() const;

    virtual std::string getPythonModule() const;

    virtual void write(OutputArchiveHandle & handle) const;

private:
    int _nonlinearDim;
    int _linearDim;
    Container _samples;
    PTR(Prior) _prior;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_SampleSet_h_INCLUDED
