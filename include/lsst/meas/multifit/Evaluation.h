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

#ifndef LSST_MEAS_MULTIFIT_Evaluation
#define LSST_MEAS_MULTIFIT_Evaluation

#include "lsst/meas/multifit/BaseEvaluator.h"

namespace lsst { namespace meas { namespace multifit {

class Evaluation : private boost::noncopyable {
public:

    Evaluation(
        BaseEvaluator::Ptr const & evaluator,
        lsst::ndarray::Array<double const,1,1> const & parameters
    );

    Evaluation(
        BaseEvaluator::Ptr const & evaluator,
        lsst::ndarray::Array<double const,1,1> const & parameters,
        BaseDistribution const & prior
    );

    Evaluation(
        BaseEvaluator::Ptr const & evaluator,
        Eigen::VectorXd const & parameters
    );

    Evaluation(
        BaseEvaluator::Ptr const & evaluator,
        Eigen::VectorXd const & parameters,
        BaseDistribution const & prior
    );

    void update(lsst::ndarray::Array<double const,1,1> const & parameters);
    void update(Eigen::VectorXd const & parameters);

    BaseEvaluator::Ptr getEvaluator() const { return _evaluator; }

    BaseDistribution::ConstPtr getPrior() const { return _prior; }
    
    lsst::ndarray::Array<double const,1,1> getParameters() const { return _parameters; }

    lsst::ndarray::Array<double const,2,2> getModelMatrix() const {
        ensure(HAS_MODEL_MATRIX);
        return _modelMatrix;
    }

    lsst::ndarray::Array<double const,3,3> getModelMatrixDerivative() const {
        ensure(HAS_MODEL_MATRIX_DERIVATIVE);
        return _modelMatrix;
    }

    lsst::ndarray::Array<double const,1,1> getCoefficients() const {
        ensure(HAS_COEFFICIENTS);
        return _coefficients;
    }

    lsst::ndarray::Array<double const,1,1> getResiduals() const {
        ensure(HAS_RESIDUALS);
        return _residuals;
    }

    lsst::ndarray::Array<double const,2,2> getResidualsJacobian() const {
        ensure(HAS_RESIDUALS_JACOBIAN);
        return _residualsJacobian;
    }

    lsst::ndarray::Array<double const,2,2> getCoefficientFisherMatrix() const {
        ensure(HAS_COEFFICIENT_FISHER_MATRIX);
        return _coefficientFisherMatrix;
    }

    lsst::ndarray::Array<double const,2,2> getCoefficientFisherFactor() const {
        ensure(HAS_COEFFICIENT_FISHER_FACTOR);
        return _coefficientFisherFactor;
    }

    double getLogPosterior() const {
        ensure(HAS_LOG_POSTERIOR);
        return _logPosterior;
    }

    double getMarginalLogPosterior() const {
        ensure(HAS_MARGINAL_LOG_POSTERIOR);
        return _marginalLogPosterior;
    }
        
private:

    class LinearSolver;

    enum StatusFlags {
        HAS_MODEL_MATRIX               = 1 <<  0,
        HAS_MODEL_MATRIX_DERIVATIVE    = 1 <<  1,
        HAS_COEFFICIENTS               = 1 <<  2,
        HAS_ACTIVE_SOLVER              = 1 <<  3,
        HAS_RESIDUALS                  = 1 <<  4,
        HAS_RESIDUALS_JACOBIAN         = 1 <<  5,
        HAS_LOG_POSTERIOR              = 1 <<  6,
        HAS_MARGINAL_LOG_POSTERIOR     = 1 <<  7
    };

    void ensure(int status) const;

    void checkSizes() const;

    int _status;
    BaseEvaluator::Ptr _evaluator;
    BaseDistribution::ConstPtr _prior;
    boost::scoped_ptr<LinearSolver> _solver;
    double _logPosterior;
    double _marginalLogPosterior;
    ndarray::Array<double,1,1> _parameters;
    mutable ndarray::Array<double,2,2> _modelMatrix;
    mutable ndarray::Array<double,3,3> _modelMatrixDerivative;
    mutable ndarray::Array<double,1,1> _coefficients;
    mutable ndarray::Array<double,2,2> _coefficientFisherMatrix;
    mutable ndarray::Array<double,2,2> _coefficientFisherFactor;
    mutable ndarray::Array<double,1,1> _residuals;
    mutable ndarray::Array<double,2,2> _residualsJacobian;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_Evaluation
