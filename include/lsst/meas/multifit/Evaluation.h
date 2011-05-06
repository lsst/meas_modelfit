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
#include "lsst/meas/multifit/GaussianDistribution.h"

namespace lsst { namespace meas { namespace multifit {

class Evaluation : private boost::noncopyable {
public:

    Evaluation(
        BaseEvaluator::Ptr const & evaluator
    );

    Evaluation(
        BaseEvaluator::Ptr const & evaluator,
        BaseDistribution const & prior
    );

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
    
    void update(
        lsst::ndarray::Array<double const,1,1> const & parameters, 
        lsst::ndarray::Array<double const,1,1> const & coefficients
    );
    void update(
        Eigen::VectorXd const & parameters,
        Eigen::VectorXd const & coefficients
    );

    void setCoefficients(lsst::ndarray::Array<double const,1,1> const & coefficients);
    void setCoefficients(Eigen::VectorXd const & coefficients);
    
    void solveCoefficients();

    BaseEvaluator::Ptr getEvaluator() const { return _evaluator; }

    BaseDistribution::ConstPtr getPrior() const { return _prior; }
    
    lsst::ndarray::Array<double const,1,1> getParameters() const { return _parameters; }

    lsst::ndarray::Array<double const,2,2> getModelMatrix() const {
        ensureModelMatrix();
        return _modelMatrix;
    }

    lsst::ndarray::Array<double const,3,3> getModelMatrixDerivative() const {
        ensureModelMatrixDerivative();
        return _modelMatrixDerivative;
    }

    lsst::ndarray::Array<double const,1,1> getCoefficients() const {
        ensureCoefficients();
        return _coefficients;
    }

    lsst::ndarray::Array<double const,1,1> getResiduals() const {
        ensureResiduals();
        return _residuals;
    }

    lsst::ndarray::Array<double const,2,2> getResidualsJacobian() const {
        ensureResidualsJacobian();
        return _residualsJacobian;
    }

    lsst::ndarray::Array<double const,2,2> getCoefficientFisherMatrix() const {
        ensureCoefficientFisherMatrix();
        return _coefficientFisherMatrix;
    }

    lsst::ndarray::Array<double const,2,2> getCoefficientFisherFactor() const {
        ensureCoefficientFisherFactor();
        return _coefficientFisherFactor;
    }

    double getLogPosterior() const {
        ensureLogPosterior();
        return _logPosterior;
    }

    ~Evaluation();

private:

#ifndef SWIG
    class LinearSolver;
    class CholeskySolver;
    // TODO: a QR solver for when we don't have a prior on the coefficients.  
    // But Eigen 2's QR solver doesn't have what we need.
#endif

    void ensureModelMatrix() const;
    void ensureModelMatrixDerivative() const;
    void ensureCoefficients() const;
    void ensureResiduals() const;
    void ensureResidualsJacobian() const;
    void ensureCoefficientFisherMatrix() const;
    void ensureCoefficientFisherFactor() const;
    void ensureLogPosterior() const;

    void initialize();
    void updateNestedPrior();

    mutable int _status;
    BaseEvaluator::Ptr _evaluator;
    BaseDistribution::ConstPtr _prior;
    GaussianDistribution::Ptr _nestedPrior;
    boost::scoped_ptr<LinearSolver> _solver;
    mutable double _logPosterior;
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
