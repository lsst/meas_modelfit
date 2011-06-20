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

#ifndef LSST_MEAS_MULTIFIT_BallCoefficientPrior
#define LSST_MEAS_MULTIFIT_BallCoefficientPrior

#include "lsst/meas/multifit/CoefficientPrior.h"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief A subclass of CoefficientPrior that sets a bound on a weighted L2 norm of the coefficient vector.
 *
 *  A prior that subclasses BallCoefficientPrior may be more sophisticated than the simple
 *  L2 constraint, but it must evaluate to zero outside the constraint.
 */
class BallCoefficientPrior : public CoefficientPrior {
public:

    typedef boost::shared_ptr<BallCoefficientPrior> Ptr;
    typedef boost::shared_ptr<BallCoefficientPrior const> ConstPtr;

    struct Constraint {
        int offset;
        int size;
        double radius;
    };

    typedef std::vector<Constraint> ConstraintVector;
    typedef ConstraintVector::const_iterator ConstraintIter;

    /// @brief Evaluate the value of the prior for a given coefficient vector.
    virtual double operator()(lsst::ndarray::Array<Pixel const,1,1> const & coefficients) const;

    /**
     *  @brief Compute the Monte Carlo integral @f$\int d\mu\,P(x|\mu,\phi)\,P(\mu|\phi) = P(x|\phi)@f$.
     *
     *  The full likelihood @f$P(x|\mu,\phi)@f$ has the form
     *  @f[
     *      \frac{1}{2}(A\mu - x)^T(A\mu - x)
     *  @f]
     *  An additional normalization factor involving the data covariance matrix must be manually included
     *  by the user.  Subclass implementations should consider the possibility that @f$A@f$ is
     *  not full rank.
     *  
     *  @param[in]  engine        Random number generator.
     *  @param[out] coefficients  (samples)x(coefficient count) array to fill with Monte Carlo samples.
     *  @param[out] weights       (samples) array of normalized weights (sum to one).
     *  @param[in]  modelMatrix   (pixel count)x(coefficient count) model matrix @f$A@f$.
     *  @param[in]  dataVector    (pixel count) array of pixel values @f$x@f$.
     *
     *  @returns a Monte Carlo estimate of @f$P(x|\phi)@f$
     */
    virtual double integrate(
        Random & engine,
        lsst::ndarray::Array<Pixel,2,2> const & coefficients,
        lsst::ndarray::Array<Pixel,1,1> const & weights,
        lsst::ndarray::Array<double const,2,2> const & modelMatrix,
        lsst::ndarray::Array<double const,1,1> const & dataVector
    ) const;

    virtual ~BallCoefficientPrior();

protected:

    explicit BallCoefficientPrior(ConstraintVector const & constraints);

    /// @brief Like operator(), but does not check constraints.  Used to implement operator().
    virtual double evaluate(lsst::ndarray::Array<Pixel const,1,1> const & coefficients) const = 0;

    /// @brief Return the volume of a hypersphere of dimension d with the given radius.
    static double computeVolume(int d, double radius);
    
    /// @brief Return the volume of a hypersphere of dimension d with unit radius.
    static double computeVolumeFactor(int d);

    /**
     *  @brief Vector of constraints.
     *
     *  Subclasses may modify the radii of the constraints but not the other members.
     */
    ConstraintVector _constraints;

private:
    class Impl;
    boost::scoped_ptr<Impl> _impl;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_BaseEvaluator
