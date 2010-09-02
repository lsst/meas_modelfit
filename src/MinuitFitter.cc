#include "lsst/meas/multifit/MinuitFitter.h"
#include <Minuit2/FCNGradientBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Eigen/Core>
#include "lsst/pex/exceptions/Runtime.h"


namespace multifit=lsst::meas::multifit;

//anonymous namespace
namespace {    

class MinuitFunction : public ROOT::Minuit2::FCNGradientBase {
public:
    MinuitFunction(multifit::ModelEvaluator::Ptr evaluator) :
        _evaluator(evaluator),
        _sigma(evaluator->computeSigmaVector()),        
        _measured(evaluator->getDataVector().getData(), evaluator->getNPixels(), 1)
    {}

    virtual ~MinuitFunction() {}
    virtual double operator()(std::vector<double> const &params) const {
        if(hasChanged(params)) { 
            _evaluator->setLinearParameters(&params[0]);
            _evaluator->setNonlinearParameters(&params[_evaluator->getLinearParameterSize()]);
        }
        ndarray::Array<multifit::Pixel const, 1, 1> modelImage(
            _evaluator->computeModelImage()
        );
        multifit::VectorMap modeled(
            _evaluator->computeModelImage().getData(), _evaluator->getNPixels()
        );
        Eigen::VectorXd residual = (_measured-modeled).cwise()/_sigma;
        return 0.5*residual.dot(residual);
    }
    virtual bool CheckGradient() const {return false;}
    virtual double Up() const {return 1.0;}
    virtual std::vector<double> Gradient(std::vector<double> const &params) const {
        if(hasChanged(params)) {
            _evaluator->setLinearParameters(&params[0]);
            _evaluator->setNonlinearParameters(&params[_evaluator->getLinearParameterSize()]);
        }
        multifit::VectorMap modeled(
            _evaluator->computeModelImage().getData(), _evaluator->getNPixels()
        );
        Eigen::VectorXd residual = (_measured-modeled).cwise()/_sigma;
        
        ndarray::Array<multifit::Pixel const, 2, 2> lpd(
            _evaluator->computeLinearParameterDerivative()
        );
        ndarray::Array<multifit::Pixel const, 2, 2> npd(
            _evaluator->computeNonlinearParameterDerivative()
        );

        multifit::MatrixMap lpdMap(
            lpd.getData(),
            _evaluator->getNPixels(),
            _evaluator->getLinearParameterSize()
        );
        multifit::MatrixMap npdMap(
            npd.getData(),
            _evaluator->getNPixels(),
            _evaluator->getNonlinearParameterSize()
        );
        
        std::vector<double> gradient(params.size());
        multifit::VectorMap gradMap(&gradient[0], params.size());
        gradMap << (-lpdMap).transpose()*residual, (-npdMap).transpose()*residual;

        return gradient;
    }
    
private:
    bool hasChanged(std::vector<double> const & params) const {
        std::vector<double>::const_iterator iNew(params.begin());
        
        double const * iOld = _evaluator->getLinearParameters().data();
        for(int i = 0; i < _evaluator->getLinearParameterSize(); ++i, ++iNew, ++iOld) {
            if(*iNew != *iOld)
                return true;
        }

        iOld = _evaluator->getNonlinearParameters().data();
        for(int i = 0; i < _evaluator->getNonlinearParameterSize(); ++i, ++iNew, ++iOld) {
            if(*iNew != *iOld)
                return true;
        }
        return false;
    }
    mutable multifit::ModelEvaluator::Ptr _evaluator;
    multifit::VectorMap _measured;
    Eigen::VectorXd _sigma;
};

}//end anonymous namespace

multifit::MinuitFitter::MinuitFitter(
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

multifit::MinuitFitter::Result multifit::MinuitFitter::apply(
    multifit::ModelEvaluator::Ptr evaluator, 
    std::vector<double> initialErrors
) const {
    ::MinuitFunction function(evaluator);
/*
    if(!function.CheckGradient()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Not computing gradient well enough according to Minuit"
        );
    }
  */      
    int nParams = evaluator->getLinearParameterSize();
    nParams += evaluator->getNonlinearParameterSize();

    if(nParams != initialErrors.size()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Number of model parameters not equal to length of error vector"
        );
    }

    std::vector<double> initialParams(nParams);
    VectorMap paramMap(&initialParams[0], nParams);
    paramMap << evaluator->getLinearParameters(), evaluator->getNonlinearParameters();
    ROOT::Minuit2::MnMigrad migrad(function, initialParams, initialErrors, _policy->getInt("strategy"));

    ROOT::Minuit2::FunctionMinimum min = migrad(
        _policy->getInt("iterationMax"), 
        _policy->getDouble("tolerance")
    );
    return Result(min);
}
