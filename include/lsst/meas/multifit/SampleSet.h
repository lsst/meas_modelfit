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

#include <list>

#include "lsst/afw/table/Catalog.h"
#include "lsst/afw/table/io/Persistable.h"
#include "lsst/meas/multifit/constants.h"
#include "lsst/meas/multifit/LogGaussian.h"
#include "lsst/meas/multifit/priors.h"
#include "lsst/meas/multifit/KernelDensityEstimator.h"

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
    Pixel marginal;    ///< Marginal nonnormalized log posterior @f$\ln m_n@f$
    Pixel proposal;    ///< Log density @f$\ln q_n@f$ of the distribution from which the samples were drawn
    double weight;     ///< Weight of this sample point; @f$m_n/q_n@f$
    Vector parameters; ///< Nonlinear parameters @f$\theta_n@f$ at this point

    /// Comparison used to sort SamplePoints by weight in order to avoid round-off error
    bool operator<(SamplePoint const & other) const { return weight < other.weight; }

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
    typedef std::list<SamplePoint> Container;
public:

    typedef Container::iterator iterator;
    typedef Container::const_iterator const_iterator;

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

    /// Return the number of nonlinear parameters
    int getNonlinearDim() const { return _nonlinearDim; }

    /// Return the number of linear amplitude parameters
    int getLinearDim() const { return _linearDim; }

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
        Eigen::VectorXd const & parameters,
        afw::geom::Point2D const & center=afw::geom::Point2D()
    ) const;

    /// @brief Return an afw::table::BaseCatalog representation of the SampleSet.
    afw::table::BaseCatalog asCatalog() const;

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

    /**
     *  @brief Add a new sample point to the SampleSet.
     *
     *  If a prior has already been applied to the SampleSet, new points may not
     *  be added (throws LogicErrorException).
     */
    void add(SamplePoint const & p);

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
    Eigen::VectorXd computeQuantiles(Eigen::VectorXd const & fractions, int parameterIndex) const;

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
    Eigen::MatrixXd computeQuantiles(Eigen::VectorXd const & fractions) const;

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
    Eigen::VectorXd computeExpectation(ExpectationFunctor const & functor, Eigen::MatrixXd * mcCov=0) const;

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
    Eigen::VectorXd computeMean(Eigen::MatrixXd * mcCov=0) const;

    /**
     *  @brief Compute the empirical covariance of the marginal distribution.
     *
     *  Using the notation of ExpectationFunctor, this is the expectation value for
     *  @f$f(\alpha,\theta)=\theta\theta^T@f$.
     *
     *  A prior must be attached before computeExpecation is called.
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
    std::string _ellipseType;
    PTR(Prior) _prior;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_SampleSet_h_INCLUDED
