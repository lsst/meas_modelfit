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

#ifndef LSST_MEAS_MULTIFIT_StudentDistribution
#define LSST_MEAS_MULTIFIT_StudentDistribution

#include "lsst/meas/multifit/SimpleDistribution.h"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief Multivariate Student distribution.
 *
 *  The density of the distribution with dimension @f$p@f$ and @f$\nu@f$ degrees of freedom is:
 *  @f[
 *     f(\mathbf{x})
 *         = \frac{\Gamma\left(\frac{\nu + p}{2}\right)}
 *                {\Gamma\left(\frac{\nu}{2}\right) \left|(\nu-2)\pi\mathbf{\Sigma}\right|^{1/2}}
 *           \left[1 + \frac{1}{\nu-2}(\mathbf{x}-\mathbf{\mu})^T \! \mathbf{\Sigma} \, 
 *                                     (\mathbf{x}-\mathbf{\mu}) \right]^{(\nu + p)/2}
 *  @f]
 *
 *  Note that the scale matrix @f$\mathbf{\Sigma}@f$ is the covariance in this definition, which
 *  is only valid for @f$\nu > 2@f$.
 */
class StudentDistribution : public SimpleDistribution {
public:

    typedef boost::shared_ptr<StudentDistribution> Ptr;
    typedef boost::shared_ptr<StudentDistribution const> ConstPtr;

    Ptr clone() const { return boost::static_pointer_cast<StudentDistribution>(_clone()); }

    virtual void draw(Random & engine, double * parameters) const;

    virtual double evaluate(double const * parameters) const;

    /// @brief Construct a Standard Normal StudentDistribution with mu=0 and sigma=I.
    explicit StudentDistribution(int dimensionality, int dof);

    /// @brief Construct a StudentDistribution with the given mu, sigma, and degrees of freedom.
    StudentDistribution(Eigen::VectorXd const & mu, Eigen::MatrixXd const & sigma, int dof);

    /// @brief Copy constructor.
    StudentDistribution(StudentDistribution const & other);

    /// @brief Assignment.
    StudentDistribution & operator=(StudentDistribution const & other);

    virtual void updateFromSamples(
        lsst::ndarray::Array<double const,2,1> const & parameters,
        lsst::ndarray::Array<double const,1,1> const & weights
    );

protected:

    virtual BaseDistribution::Ptr _clone() const;

    virtual void invalidate() { _cached.reset(); }

private:

    void ensureCached() const;

    struct Cached {
        double normalization;
        Eigen::MatrixXd factor; // lower-triangular Cholesky factor of sigma.
    };

    int _dof;
    mutable boost::shared_ptr<Cached> _cached;
    mutable Eigen::VectorXd _workspace;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_StudentDistribution
