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

#ifndef LSST_MEAS_MULTIFIT_BaseDistribution
#define LSST_MEAS_MULTIFIT_BaseDistribution

#include "lsst/ndarray.h"
#include "lsst/pex/exceptions.h"
#include "lsst/meas/multifit/constants.h"
#include <Eigen/Core>

namespace lsst { namespace meas { namespace multifit {

class BaseDistribution {
public:

    typedef boost::shared_ptr<BaseDistribution> Ptr;
    typedef boost::shared_ptr<BaseDistribution const> ConstPtr;

    typedef int DependencyFlags;

    /**
     *  @brief Deep-copy the distribution.
     *
     *  Subclasses should shadow this member function with another that static-downcasts the
     *  protected implementation.
     */
    Ptr clone() const { return _clone(); }

    /// @brief Return the dimensionality of the distribution.
    int getDimensionality() const { return _dimensionality; }

    ///@{
    /// @brief Draw a parameter vector from the distribution.
    virtual void draw(Random & engine, double * parameters) const = 0;
    void draw(Random & engine, lsst::ndarray::Array<double,1,1> const & parameters) const {
        draw(engine, parameters.begin());
    }
    void draw(Random & engine, Eigen::VectorXd & parameters) const {
        draw(engine, parameters.data());
    }
    ///@}

    ///@{
    ///  @brief Evaluate the distribution at the given parameters.
    virtual double evaluate(double const * parameters) const = 0;
    double evaluate(lsst::ndarray::Array<double,1,1> const & parameters) const {
        return evaluate(parameters.begin());
    }
    double evaluate(Eigen::VectorXd const & parameters) const {
        return evaluate(parameters.data());
    }
    ///@}

    /// @brief Return the dimensionality of the nested distribution.
    virtual int getNestedDimensionality() const { return 0; }

    /**
     *  @brief Return a set of bitflags that defined how the nested conditional distribution
     *         depends on the parameters of the outer distribution.
     *
     *  A return value of zero indicates that the nested distribution is completely independent,
     *  making calls to updateNested() completely unnecessary.  The meaning of any nonzero
     *  value is specific to the type of the nested dependency, which may have additional methods
     *  to interpret the flags.
     */
    virtual DependencyFlags getNestedDependency() const { return 0; }

    ///@{
    /**
     *  @brief Extract a nested distribution at the given parameters.
     *
     *  Nesting allows a distribution to model a joint distribution P(x,y) as the product
     *  of its own marginal distribution P(x) and the nested conditional distribution P(y|x).
     *
     *  Most distributions are not nested and return an empty pointer.  Subclasses of BaseDistribution
     *  that do support nesting should shadow these functions with static downcasting versions that
     *  return the appropriate nested distribution subclass, but still make use of the virtual
     *  protected implementation.
     */
    Ptr evaluateNested(double const * parameters) const {
        return _evaluateNested(parameters);
    }
    Ptr evaluateNested(lsst::ndarray::Array<double,1,1> const & parameters) const {
        return _evaluateNested(parameters.getData());
    }
    Ptr evaluateNested(Eigen::VectorXd const & parameters) const {
        return _evaluateNested(parameters.data());
    }
    ///@}

    ///@{
    /**
     *  @brief Update a nested distribution to a new set of parameters.
     *
     *  Only a distribution returned by evaluateNested() should be passed to updateNested(); passing
     *  in another distribution may result in undefined behavior.
     *
     *  The default implementation throws an exception.
     */
    void updateNested(BaseDistribution & nested, double const * parameters) const {
        _updateNested(nested, parameters);
    }
    void updateNested(BaseDistribution & nested, lsst::ndarray::Array<double,1,1> const & parameters) const {
        _updateNested(nested, parameters.getData());
    }
    void updateNested(BaseDistribution & nested, Eigen::VectorXd const & parameters) const {
        _updateNested(nested, parameters.data());
    }
    ///@}

    /// @brief Compute the mean of the distribution.
    virtual Eigen::VectorXd computeMean() const = 0;

    /// @brief Return the covariance matrix of the distribution.
    virtual Eigen::MatrixXd computeCovariance() const = 0;

    /**
     *  @brief Modify the distribution to match a set of importance or MCMC samples.
     *
     *  This will generally be used by adaptive importance sampling methods, and most
     *  operations will match moments or minimize the Kullback-Leibler divergence.
     *
     *  The default implementation throws an exception.
     *
     *  The number of rows of the parameters matrix should match the number of elements in the
     *  weights vector, while the number of columns matches the dimensionality of the distribution.
     */
    virtual void updateFromSamples(
        lsst::ndarray::Array<double const,2,1> const & parameters,
        lsst::ndarray::Array<double const,1,1> const & weights
    ) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "updateFromSamples is not implemented for this distribution."
        );
    }

protected:

    virtual Ptr _clone() const = 0;

    virtual Ptr _evaluateNested(double const * parameters) const { return Ptr(); }

    virtual void _updateNested(BaseDistribution & nested, double const * parameters) const {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "updateNested is not implemented for this distribution."
        );
    }

    BaseDistribution(int dimensionality) : _dimensionality(dimensionality) {}

    BaseDistribution(BaseDistribution const & other) : _dimensionality(other._dimensionality) {}

    void operator=(BaseDistribution const & other) { _dimensionality = other._dimensionality; }

    int _dimensionality;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_BaseDistribution
