// -*- lsst-c++ -*-
/**
 * @file 
 * Implementation of a SingleLinearParameterFitter
 */
#include <cfloat>
#include <cmath>
#include <iostream>
#include <Eigen/Array>
#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/Cholesky>

#include "lsst/meas/multifit/SingleLinearParameterFitter.h"
#include "lsst/meas/multifit/core.h"
#include "lsst/pex/exceptions/Runtime.h"

namespace multifit = lsst::meas::multifit;

/**
 * Configure a SingleLinearParameterFitter with a policy
 *
 * Some default value may be provided by the dictionary 
 */
multifit::SingleLinearParameterFitter::SingleLinearParameterFitter(
    lsst::pex::policy::Policy::Ptr const & policy
) : _policy(policy) {
    if(!_policy)
        _policy.reset(new lsst::pex::policy::Policy());

    //load default policy
    lsst::pex::policy::Policy::Ptr defaults(
        lsst::pex::policy::Policy::createPolicy(*getDefaultPolicySource())
    );
    //merge in default values
    if(defaults->canValidate()){
        _policy->mergeDefaults(*defaults->getDictionary());
    }
    else {
        _policy->mergeDefaults(*defaults);
    }


    if(!_policy->exists("terminationType")) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Invalid configuration policy - missing value \"terminationType\""
        );
    } 
    
    _terminationType = 0;
    std::vector<std::string> conditions(_policy->getStringArray("terminationType"));
    for(std::vector<std::string>::iterator i(conditions.begin()), end(conditions.end());
        i != end; ++i
    ) {
        if((*i) == "dChisq")
            _terminationType |= DCHISQ;
        else if((*i) == "iteration")
            _terminationType |= ITERATION;
        else if((*i) == "step")
            _terminationType |= STEP; 
    }
    
    if(_terminationType & ITERATION) {
        if(!_policy->exists("iterationMax")) {
            throw LSST_EXCEPT(
                lsst::pex::exceptions::InvalidParameterException,
                "Invalid configuration policy - missing value \"iterationMax\""
            );
        } else {
            _iterationMax = _policy->getInt("iterationMax");   
        }
    }
    if(_terminationType & DCHISQ) {
        if(!_policy->exists("dChisqThreshold")) {
            throw LSST_EXCEPT(
                lsst::pex::exceptions::InvalidParameterException,
                "Invalid configuration policy - missing value \"dChisqThreshold\""
            );
        } else {
            _dChisqThreshold = std::abs(_policy->getDouble("dChisqThreshold"));   
        }
    }
    if(_terminationType & STEP) {
        if(!_policy->exists("stepThreshold")) {
            throw LSST_EXCEPT(
                lsst::pex::exceptions::InvalidParameterException,
                "Invalid configuration policy - missing value \"stepThreshold\""
            );
        } else {
            _stepThreshold = std::abs(_policy->getDouble("stepThreshold"));   
        }
    }
}    

/**
 * Fit a model by applying this to a ModelEvaluator
 *
 * Use this configured fitter to compute the best-fit model for the given
 * ModelEvaluator. The ModelEvluator must be properly initialized with the set
 * of exposures to fit a model on. 
 *
 * The fitting loop will terminated when at least one of the termination 
 * conditions specified to this fitter upon construction is met. 
 *
 * @param evaluator must be properly initialized externally, by setting
 *      exposure list upon its construction, or by calling 
 *      ModelEvaluator::setExposureList
 * @return SimpleFitResult which contains information about the status of the
 *      model fit upon termination.
 *
 * @sa lsst::meas::multifit::ModelEvaluator
 * @sa lsst::meas::multifit::SimpleFitResult
 */
multifit::SingleLinearParameterFitter::Result::Ptr multifit::SingleLinearParameterFitter::apply(
    ModelEvaluator & evaluator
) const {
    if(evaluator.getLinearParameterSize() != 1) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "SingleLinearParameterFitter can only be applied to evaluators with 1 linear Parameter"
        );
    }
    else if (evaluator.getNPixels() <= 0) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "ModelEvaluator has no associated exposures."
        );
    }
    Result::Ptr result = boost::make_shared<Result>();
    result->model = evaluator.getModel();    

    VectorMap data(
        evaluator.getDataVector().getData(),
        evaluator.getNPixels()
    );    
    Eigen::VectorXd sigma(evaluator.computeSigmaVector());

    Eigen::VectorXd observed = (data.cwise() / sigma);
    Eigen::VectorXd newLinearDNonlinear;
    Eigen::MatrixXd jacobian;


    int nIterations = 0;
    double chisq=DBL_MAX, dChisq=DBL_MAX; 
    Eigen::VectorXd step;
    bool done = false;
    while( !done) {  
        ndarray::Array<Pixel const, 2, 1> dLinearArray(evaluator.computeLinearParameterDerivative());

        VectorMap dLinear = VectorMap(
            evaluator.computeLinearParameterDerivative().getData(),
            evaluator.getNPixels()
        );
        dLinear.cwise() /= sigma;

        MatrixMap dNonlinear = MatrixMap(
            evaluator.computeNonlinearParameterDerivative().getData(),
            evaluator.getNPixels(),
            evaluator.getNonlinearParameterSize()
        );
        for(int i = 0; i < dNonlinear.cols(); ++i) {
            dNonlinear.col(i).cwise() /= sigma;
        }

        double normDLinear = dLinear.squaredNorm();        
        
        //the new linear parameter as a function of the nonlinear parameters
        double newLinear = dLinear.dot(observed) / normDLinear;
        //compute the residual as a function of the nonlinear parameters
        Eigen::VectorXd residual = observed - dLinear*newLinear;

        //compute the chisq
        if (nIterations > 0) {
            dChisq = chisq;
        }
        chisq = (residual.dot(residual))/2;
        dChisq -= chisq;
        

        //compute derivative of new linear parameter w.r.t nonlinear parameters
        //this is a matrix with dimensions (nNonlinear, 1)
        Eigen::VectorXd dNewLinear = dNonlinear.transpose()*(residual - dLinear*newLinear);
        dNewLinear /= (normDLinear*newLinear);
    
        //compute the jacobian of partial derivatives of the model w.r.t
        //nonlinear parameters
        //this is a matrix with dimensions (pixels, nNonlinear)
        Eigen::MatrixXd jacobian = -dNonlinear - dLinear*dNewLinear.transpose();

        //compute the step to take on nonlinear parameters:
        //this is a matrix with dimensions (nNonlinear, 1)
        step = jacobian.transpose()*residual;

        (jacobian.transpose()*jacobian).llt().solveInPlace(step);

        Eigen::VectorXd newNonlinear = evaluator.getNonlinearParameters() + step;
        evaluator.setLinearParameters(&newLinear);
        evaluator.setNonlinearParameters(newNonlinear.data());
        ++nIterations;
        
        if (_terminationType & ITERATION && nIterations >= _iterationMax) {
            done = true;
            result->convergenceFlags |= Result::MAX_ITERATION_REACHED;
        }
        if( (_terminationType & DCHISQ) && (nIterations > 1) && (std::abs(dChisq) < _dChisqThreshold) ) {
            done = true; 
            result->convergenceFlags |= Result::DCHISQ_THRESHOLD_REACHED;
            result->convergenceFlags |= Result::CONVERGED;
        }
        if( (_terminationType & STEP) && (nIterations > 1) && (step.norm() < _stepThreshold) ) {
            done = true;
            result->convergenceFlags |= Result::STEP_THRESHOLD_REACHED;
            result->convergenceFlags |= Result::CONVERGED; 
        }
    };

    result->chisq = chisq;
    result->dChisq = dChisq;
    result->sdqaMetrics->set("nIterations", nIterations);
    result->sdqaMetrics->set("finalNonlinearStepNorm", step.norm());

    return result;
}

