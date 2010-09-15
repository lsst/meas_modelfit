#include "lsst/meas/multifit/MinuitFitter.h"
#include <Minuit2/FCNGradientBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/FunctionGradient.h>
#include <Minuit2/LAVector.h>
#include <Minuit2/InitialGradientCalculator.h>
#include <Minuit2/Numerical2PGradientCalculator.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MinimumParameters.h>
#include <Minuit2/MnFcn.h>
#include <Eigen/Core>
#include "lsst/pex/exceptions/Runtime.h"
#include <sstream>
#include "lsst/pex/logging/Trace.h"

namespace multifit=lsst::meas::multifit;
namespace pexLog = lsst::pex::logging;

//anonymous namespace
namespace {   

class ChisqFunction {
public:
    typedef Eigen::Matrix<multifit::Pixel, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Matrix<multifit::Pixel, Eigen::Dynamic, Eigen::Dynamic> Matrix;

    ChisqFunction(
        multifit::ModelEvaluator::Ptr const & evaluator,
        Vector const & priorMean,
        Vector const & priorFisherDiag
    ) 
      : _dirty(true), 
        _evaluator(evaluator),
        _measured(evaluator->getWeightedData()),
        _priorMean(priorMean),
        _priorFisherDiag(priorFisherDiag)
    {}

    
    double computeValue(std::vector<double> const &params) {
        checkParams(params);
        if(_dirty) {
            _evaluator->setLinearParameters(&params[0]);
            _evaluator->setNonlinearParameters(&params[_evaluator->getLinearParameterSize()]);
        }

        Vector modeled = _evaluator->computeModelImage();
        Vector residual = (_measured-modeled);
        Vector priorDelta(params.size());
        std::copy(params.begin(), params.end(), priorDelta.data());
        priorDelta -= _priorMean;
        double priorTerm = 0.5*priorDelta.dot(_priorFisherDiag.cwise() * priorDelta);
        return 0.5*residual.dot(residual) + priorTerm; 
    }

    std::vector<double> computeGradient(std::vector<double> const &params) {
        checkParams(params);
        if(_dirty) {
            _evaluator->setLinearParameters(&params[0]);
            _evaluator->setNonlinearParameters(&params[_evaluator->getLinearParameterSize()]);
        }

        Vector modeled = _evaluator->computeModelImage();
        Vector residual = (_measured-modeled);
        Matrix lpd = _evaluator->computeLinearParameterDerivative();
        Matrix npd = _evaluator->computeNonlinearParameterDerivative();

        Eigen::VectorXd priorDelta(params.size());
        std::copy(params.begin(), params.end(), priorDelta.data());
        priorDelta -= _priorMean;

        std::vector<double> gradient(params.size());
        multifit::VectorMap gradMap(&gradient[0], params.size());
        gradMap << (-lpd).transpose()*residual, (-npd).transpose()*residual;
        gradMap += _priorFisherDiag.cwise() * priorDelta;
        return gradient; 
    }
private:
    bool _dirty;
    multifit::ModelEvaluator::Ptr _evaluator;
    Vector _measured;
    Vector _priorMean;
    Vector _priorFisherDiag;

    void checkParams(std::vector<double> const & params) {
        if(_dirty)
            return;

        std::vector<double>::const_iterator iNew(params.begin());
	       
        double const * iOld = _evaluator->getLinearParameters().data();
        for(int i = 0; i < _evaluator->getLinearParameterSize(); ++i, ++iNew, ++iOld) {
            if(*iNew != *iOld)
                _dirty = true;
                return;
        }

        iOld = _evaluator->getNonlinearParameters().data();
        for(int i = 0; i < _evaluator->getNonlinearParameterSize(); ++i, ++iNew, ++iOld) {
            if(*iNew != *iOld)
                _dirty = true;
                return;
        }

        _dirty = false;
    }
};


class Function : public ROOT::Minuit2::FCNBase {
public:
    Function(
        multifit::ModelEvaluator::Ptr evaluator, 
        ChisqFunction::Vector const & priorMean,
        ChisqFunction::Vector const & priorFisherDiag
    ) :
        _chisqFunction(evaluator, priorMean, priorFisherDiag)
    {}

    virtual ~Function() {}
    virtual double operator()(std::vector<double> const &params) const {
        return _chisqFunction.computeValue(params);
    }
    virtual double Up() const {return 1.0;}
private:
    mutable ChisqFunction _chisqFunction;
};


class GradientFunction : public ROOT::Minuit2::FCNGradientBase {
public:
    GradientFunction(
        multifit::ModelEvaluator::Ptr evaluator,
        ChisqFunction::Vector const & priorMean,
        ChisqFunction::Vector const & priorFisherDiag,
        bool checkGradient=false
    ) :
        _chisqFunction(evaluator, priorMean, priorFisherDiag),
        _checkGradient(checkGradient)
    {}

    virtual ~GradientFunction() {}
    virtual double operator()(std::vector<double> const &params) const {
        return _chisqFunction.computeValue(params);
    }
    virtual bool CheckGradient() const {return _checkGradient;}
    virtual double Up() const {return 1.0;}
    virtual std::vector<double> Gradient(std::vector<double> const &params) const {
        return _chisqFunction.computeGradient(params);
    }
private:
    mutable ChisqFunction _chisqFunction;
    bool _checkGradient;
};

}//end anonymous namespace


/******************************************************************************
 * Analytic Fitter
 *****************************************************************************/
multifit::MinuitAnalyticFitter::MinuitAnalyticFitter(
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
}

multifit::MinuitFitterResult multifit::MinuitAnalyticFitter::apply(
    multifit::ModelEvaluator::Ptr evaluator, std::vector<double> const & initialErrors
) const {
    int nParams = evaluator->getLinearParameterSize() + evaluator->getNonlinearParameterSize();
    std::vector<double> priorMean(nParams, 0.0);
    std::vector<double> priorFisherDiag(nParams, 0.0);
    return apply(evaluator, initialErrors, priorMean, priorFisherDiag);
}

multifit::MinuitFitterResult multifit::MinuitAnalyticFitter::apply(
    multifit::ModelEvaluator::Ptr evaluator, 
    std::vector<double> const & initialErrors,
    std::vector<double> const & priorMean,
    std::vector<double> const & priorFisherDiag
) const {
    
    int nParams = evaluator->getLinearParameterSize();
    nParams += evaluator->getNonlinearParameterSize();

    if(nParams != initialErrors.size()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Number of model parameters not equal to length of error vector"
        );
    }

    bool checkGradient = _policy->getBool("checkGradient");

    ChisqFunction::Vector priorMeanEig(priorMean.size());
    std::copy(priorMean.begin(), priorMean.end(), priorMeanEig.data());

    ChisqFunction::Vector priorFisherDiagEig(priorFisherDiag.size());
    std::copy(priorFisherDiag.begin(), priorFisherDiag.end(), priorFisherDiagEig.data());

    if(nParams != priorMean.size() || nParams != priorFisherDiag.size()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Number of model parameters not equal to length of prior mean and/or Fisher diagonal vector"
        );
    }

    ::GradientFunction function(evaluator, priorMeanEig, priorFisherDiagEig);

    std::vector<double> initialParams(nParams);
    VectorMap paramMap(&initialParams[0], nParams);
    paramMap << evaluator->getLinearParameters(), evaluator->getNonlinearParameters();

    ROOT::Minuit2::MnMigrad migrad(function, initialParams, initialErrors, _policy->getInt("strategy"));

    if(checkGradient) {
        ROOT::Minuit2::MnFcn mnFcn(function);
        ROOT::Minuit2::LAVector paramVec(nParams);
        for(int i =0; i< nParams; ++i) {
            paramVec[i] =  initialParams[i];
        }   

        ROOT::Minuit2::MinimumParameters minParam(paramVec, function(initialParams));
        ROOT::Minuit2::FunctionGradient fncGrad = ROOT::Minuit2::Numerical2PGradientCalculator(
                mnFcn, 
                migrad.State().Trafo(), 
                migrad.Strategy()
            )(minParam);

        ROOT::Minuit2::LAVector numericGrad = fncGrad.Grad();        
        std::vector<double> analyticGrad = function.Gradient(initialParams);
        std::ostringstream numericStr, analyticStr, diffStr;

        numericStr << "numeric gradient: <" << numericGrad[0];
        analyticStr << "analytic gradient: <" << analyticGrad[0];
        diffStr << "ratio: <" << numericGrad[0]/analyticGrad[0];

        for( int i =1; i < nParams; ++i) {
            numericStr << ", " << numericGrad[i];
            analyticStr << ", " << analyticGrad[i];
            diffStr << ", " << numericGrad[i]/analyticGrad[i];
        }
        numericStr << ">";
        analyticStr << ">";
        diffStr << ">";

        std::cerr << numericStr.str() <<std::endl;
        std::cerr << analyticStr.str() <<std::endl;
        std::cerr << diffStr.str() <<std::endl;
    }
   
    ROOT::Minuit2::FunctionMinimum min = migrad(
        _policy->getInt("iterationMax"), 
        _policy->getDouble("tolerance")
    );
    return Result(min, evaluator->getModel());
}

/******************************************************************************
 * Numeric Fitter
 *****************************************************************************/
multifit::MinuitNumericFitter::MinuitNumericFitter(
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
}

multifit::MinuitFitterResult multifit::MinuitNumericFitter::apply(
    multifit::ModelEvaluator::Ptr evaluator, std::vector<double> const & initialErrors
) const {
    int nParams = evaluator->getLinearParameterSize() + evaluator->getNonlinearParameterSize();
    std::vector<double> priorMean(nParams, 0.0);
    std::vector<double> priorFisherDiag(nParams, 0.0);
    return apply(evaluator, initialErrors, priorMean, priorFisherDiag);
}

multifit::MinuitFitterResult multifit::MinuitNumericFitter::apply(
    multifit::ModelEvaluator::Ptr evaluator, 
    std::vector<double> const & initialErrors,
    std::vector<double> const & priorMean,
    std::vector<double> const & priorFisherDiag
) const {
    
    int nParams = evaluator->getLinearParameterSize();
    nParams += evaluator->getNonlinearParameterSize();

    if(nParams != initialErrors.size()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Number of model parameters not equal to length of error vector"
        );
    }

    Eigen::VectorXd priorMeanEig(priorMean.size());
    std::copy(priorMean.begin(), priorMean.end(), priorMeanEig.data());
    Eigen::VectorXd priorFisherDiagEig(priorFisherDiag.size());
    std::copy(priorFisherDiag.begin(), priorFisherDiag.end(), priorFisherDiagEig.data());

    if(nParams != priorMean.size() || nParams != priorFisherDiag.size()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Number of model parameters not equal to length of prior mean and/or Fisher diagonal vector"
        );
    }

    ::Function function(evaluator, priorMeanEig, priorFisherDiagEig);
    
    std::vector<double> initialParams(nParams);
    VectorMap paramMap(&initialParams[0], nParams);
    paramMap << evaluator->getLinearParameters(), evaluator->getNonlinearParameters();

    ROOT::Minuit2::MnMigrad migrad(function, initialParams, initialErrors, _policy->getInt("strategy"));

    ROOT::Minuit2::FunctionMinimum min = migrad(
        _policy->getInt("iterationMax"), 
        _policy->getDouble("tolerance")
    );

    return Result(min, evaluator->getModel());
}
