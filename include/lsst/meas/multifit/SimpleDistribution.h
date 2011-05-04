// -*- LSST-C++ -*-
/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010, 2011 LSST Corporation.
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

#ifndef LSST_MEAS_MULTIFIT_SimpleDistribution
#define LSST_MEAS_MULTIFIT_SimpleDistribution

#include "lsst/meas/multifit/BaseDistribution.h"

namespace lsst { namespace meas { namespace multifit {

class SimpleInterpreter;

/**
 *  @brief Intermediate base class for istributions defined by a location (mu) vector
 *         and symmetric positive-definite covariance (sigma) matrix.
 */
class SimpleDistribution : public BaseDistribution {
public:

    typedef boost::shared_ptr<SimpleDistribution> Ptr;
    typedef boost::shared_ptr<SimpleDistribution const> ConstPtr;

    Ptr clone() const { return boost::static_pointer_cast<SimpleDistribution>(_clone()); }

    /// @brief Compute the mean of the distribution.
    virtual Eigen::VectorXd computeMean() const { return _mu; }

    /// @brief Return the covariance matrix of the distribution.
    virtual Eigen::MatrixXd computeCovariance() const { return _sigma; }

    /// @brief Return the location parameters of the distribution as a vector.
    Eigen::VectorXd const & getMu() const { return _mu; }

    /// @brief Set the entire vector of location parameters at once.
    void setMu(Eigen::VectorXd const & mu) { checkShape(mu, _sigma); _mu = mu; }

    /// @brief Return the covariance of the distribution as a symmetric positive-definite matrix.
    Eigen::MatrixXd const & getSigma() const { return _sigma; }

    /// @brief Set the entire covariance matrix at once.
    void setSigma(Eigen::MatrixXd const & sigma) { checkShape(_mu, sigma); invalidate(); _sigma = sigma; }

protected:

    /// @brief Construct a SimpleDistribution with no grid with mu=0 and sigma=I.
    explicit SimpleDistribution(int dimensionality) : BaseDistribution(dimensionality) {}

    /// @brief Copy constructor.
    SimpleDistribution(SimpleDistribution const & other) : BaseDistribution(other) {}

    /// @brief Assignment.  Protected to prevent slicing.
    void operator=(SimpleDistribution const & other) {
        BaseDistribution::operator=(other);
    }

    virtual void invalidate() = 0;

    void checkShape(Eigen::VectorXd const & mu, Eigen::MatrixXd const & sigma) {
        if (sigma.rows() != mu.size() || sigma.cols() != mu.size()) {
            throw LSST_EXCEPT(
                lsst::pex::exceptions::InvalidParameterException,
                (boost::format("Sigma matrix shape (%d, %d) does not match mu vector size (%d).")
                 % sigma.rows() % sigma.cols() % mu.size()).str()
            );
        }
    }

    Eigen::VectorXd _mu;
    Eigen::MatrixXd _sigma;

private:

    friend class SimpleInterpreter;

};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_SimpleDistribution
