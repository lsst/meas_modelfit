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

#ifndef LSST_MEAS_MULTIFIT_GaussianDistribution
#define LSST_MEAS_MULTIFIT_GaussianDistribution

#include "lsst/meas/multifit/SimpleDistribution.h"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief A Gaussian or Normal distribution.
 */
class GaussianDistribution : public SimpleDistribution {
public:

    typedef boost::shared_ptr<GaussianDistribution> Ptr;
    typedef boost::shared_ptr<GaussianDistribution const> ConstPtr;

    Ptr clone() const { return boost::static_pointer_cast<GaussianDistribution>(_clone()); }

    virtual void draw(RandomEngine & engine, double * parameters) const;

    virtual double evaluate(double const * parameters) const;

    /// @brief Construct a Standard Normal GaussianDistribution with mu=0 and sigma=I.
    explicit GaussianDistribution(int size);

    /// @brief Construct a GaussianDistribution with the given mu and sigma.
    GaussianDistribution(Eigen::VectorXd const & mu, Eigen::MatrixXd const & sigma);

    /// @brief Copy constructor.
    GaussianDistribution(GaussianDistribution const & other);

    /// @brief Assignment.
    GaussianDistribution & operator=(GaussianDistribution const & other);

    virtual void updateFromSamples(
        Eigen::MatrixXd const & parameters,
        Eigen::VectorXd const & weights
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

    mutable boost::shared_ptr<Cached> _cached;
    mutable Eigen::VectorXd _workspace;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_GaussianDistribution
