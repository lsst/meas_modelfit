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

    ndarray::Array<Pixel const, 2, 2> dLinearArray, dNonlinearArray;    
    Eigen::VectorXd linearParam = evaluator.getLinearParameters();
    Eigen::VectorXd nonlinearParam = evaluator.getNonlinearParameters();

    double chisq=DBL_MAX, dChisq=DBL_MAX; 
    Eigen::VectorXd step;
    bool endOnIteration = ((_terminationType & ITERATION) != 0);
    bool done = false;
    int nIterations = 0;
    for(; ((nIterations < _iterationMax) && endOnIteration) || !done; ++nIterations) {
        try { 
            ndarray::shallow(dLinearArray) = 
                evaluator.computeLinearParameterDerivative();
        } catch(multifit::ParameterRangeException & e) {
            multifit::ParameterMap const & paramMap(e.getOutOfRangeParameterMap());
            multifit::ParameterMap::const_iterator i(paramMap.begin());
            linearParam[0] = (linearParam[0] + i->second)/2.0;
            evaluator.setLinearParameters(linearParam);
            ndarray::shallow(dLinearArray) = 
                evaluator.computeLinearParameterDerivative();
        }       

        VectorMap dLinear = VectorMap(
            dLinearArray.getData(),
            evaluator.getNPixels()
        );
        dLinear.cwise() /= sigma;


        try {    
            ndarray::shallow(dNonlinearArray) =
                evaluator.computeNonlinearParameterDerivative();
        } catch(multifit::ParameterRangeException & e) {
            multifit::ParameterMap const & paramMap(e.getOutOfRangeParameterMap());
            for(ParameterMap::const_iterator i(paramMap.begin()); i != paramMap.end(); ++i) {
                nonlinearParam[i->first] = (nonlinearParam[i->first] + i->second)/2.0;
            }
            evaluator.setNonlinearParameters(nonlinearParam);
            ndarray::shallow(dNonlinearArray) =
                evaluator.computeNonlinearParameterDerivative();
        }
        MatrixMap dNonlinear = MatrixMap(
            dNonlinearArray.getData(),
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
        if( (_terminationType & DCHISQ) && (nIterations > 0) && 
            (std::abs(dChisq) < _dChisqThreshold) 
        ) {
                done = true; 
                result->convergenceFlags |= Result::DCHISQ_THRESHOLD_REACHED;
                result->convergenceFlags |= Result::CONVERGED;
        } 

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
        if( (_terminationType & STEP) && (step.norm() < _stepThreshold) ) {
            done = true;
            result->convergenceFlags |= Result::STEP_THRESHOLD_REACHED;
            result->convergenceFlags |= Result::CONVERGED; 
        }

        (jacobian.transpose()*jacobian).llt().solveInPlace(step);


        linearParam = evaluator.getLinearParameters();
        try {
            evaluator.setLinearParameters(&newLinear);
        } catch(multifit::ParameterRangeException & e) {
            multifit::ParameterMap const & paramMap(e.getOutOfRangeParameterMap());
            multifit::ParameterMap::const_iterator i(paramMap.begin());
            newLinear = (linearParam[0] + i->second)/2.0;
            evaluator.setLinearParameters(&newLinear);
        }

        nonlinearParam = evaluator.getNonlinearParameters();
        Eigen::VectorXd newNonlinear = nonlinearParam + step;
        try {
            evaluator.setNonlinearParameters(newNonlinear);
        } catch(multifit::ParameterRangeException & e) {
            multifit::ParameterMap const & paramMap(e.getOutOfRangeParameterMap());
            for(ParameterMap::const_iterator i(paramMap.begin()); i != paramMap.end(); ++i) {
                newNonlinear[i->first] = (nonlinearParam[i->first] + i->second)/2.0;
            }
            evaluator.setNonlinearParameters(newNonlinear);
        }
    };
    
    if (nIterations >= _iterationMax) {
        result->convergenceFlags |= Result::MAX_ITERATION_REACHED;
    }
    
    result->chisq = chisq;
    result->dChisq = dChisq;
    result->sdqaMetrics->set("nIterations", nIterations);
    result->sdqaMetrics->set("finalNonlinearStepNorm", step.norm());

    return result;
}

