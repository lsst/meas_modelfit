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
#ifndef LSST_MEAS_MULTIFIT_BaseObjective
#define LSST_MEAS_MULTIFIT_BaseObjective

#include "lsst/meas/multifit/constants.h"
#include <Eigen/Core>

namespace lsst { namespace meas { namespace multifit {

class BaseObjective : private boost::noncopyable {
public:

    typedef boost::shared_ptr<Objective> Ptr;
    typedef boost::shared_ptr<Objective const> Ptr;

    /// @brief Return the number of parameters defined by the objective.
    int getParameterSize() const { return _parameters.size(); }

    /// @brief Return the current parameters.
    Eigen::VectorXd const & getParameters() const { return _parameters; }

    /// @brief Return the step that sets the test parameters.
    Eigen::VectorXd const & getStep() const { return _step; }

    /// @brief Return the value of the objective function at the current parameters. 
    double getValue() const { return _value; }

    /// @brief Return the value of the objective function at the test parameters.
    double getTestValue() const { return _testValue; }

    /// @brief Return the gradient of the objective at the current (non-test) parameters.
    virtual Eigen::VectorXd const & getGradient() const {
        throw DerivativeNotImplementedError("Gradient not implemented by objective.");
    }

    /// @brief Return the Hessian of the objective at the current (non-test) parameters.
    virtual Eigen::MatrixXd const & getHessian() const {
        throw DerivativeNotImplementedError("Hessian not implemented by objective.");
    }

    /// @brief Return the residuals vector at the current parameters.
    Eigen::VectorXd const & getResiduals() const { return _residuals; }

    /// @brief Return the residuals vector at the test parameters.
    Eigen::VectorXd const & getTestResiduals() const { return _testResiduals; }

    /// @brief Return the Jacobian of the residuals at the current (non-test) parameters.
    virtual Eigen::MatrixXd const & getJacobian() const {
        throw DerivativeNotImplementedError("Gradient not implemented by objective.");
    }

    /**
     *  @brief Set parameters = parameters + step and step = 0.  Update derivatives
     *         accordingly.
     */
    virtual void acceptStep() = 0;

    /**
     *  @brief Set the step to the given vector and return the resulting test parameter
     *         objective value.
     */
    virtual double test(Eigen::VectorXd const & step) = 0;

    virtual ~Objective() {}

protected:

    double _value;
    double _testValue;
    Eigen::VectorXd _residuals;
    Eigen::VectorXd _testResiduals;
    Eigen::VectorXd _parameters;
    Eigen::VectorXd _step;
};

} // namespace modeling

#endif // !LSST_MEAS_MULTIFIT_BaseObjective
