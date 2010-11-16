#include "lsst/meas/multifit/MinuitFitter.h"
#include <Minuit2/FCNGradientBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/FunctionGradient.h>
#include <Minuit2/LAVector.h>
#include <Minuit2/InitialGradientCalculator.h>
#include <Minuit2/Numerical2PGradientCalculator.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnSimplex.h>
#include <Minuit2/MinimumParameters.h>
#include <Minuit2/MnFcn.h>
#include <Eigen/Core>
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/pex/logging/Trace.h"
#include <limits>
#include <iostream>
#include <boost/scoped_ptr.hpp>

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
        Vector const & priorFisherDiag,
        bool doSnapshots,
        std::string const & snapshotFormat
    ) 
      : _dirty(true), 
        _doSnapshots(doSnapshots),
        _nUpdates(0),
        _snapshotFormat(snapshotFormat),
        _evaluator(evaluator),
        _measured(evaluator->getWeightedData()),
        _priorMean(priorMean),
        _priorFisherDiag(priorFisherDiag)
    {}
    
    double computeValue(std::vector<double> const &params) {
        std::cerr << "Computing value: ";
        std::copy(params.begin(), params.end(), std::ostream_iterator<double>(std::cerr, " "));
        std::cerr << std::endl;
        setParams(params);
        Vector modeled = _evaluator->computeModelImage();
        Vector residual = (_measured-modeled);
        Vector priorDelta(params.size());
        std::copy(params.begin(), params.end(), priorDelta.data());
        priorDelta -= _priorMean;
        double priorTerm = 0.5*priorDelta.dot(_priorFisherDiag.cwise() * priorDelta);
        double value = 0.5*residual.dot(residual) + priorTerm; 
        return value;
    }

    std::vector<double> computeGradient(std::vector<double> const &params) {
        std::cerr << "Computing gradient: ";
        std::copy(params.begin(), params.end(), std::ostream_iterator<double>(std::cerr, " "));
        std::cerr << std::endl;
        setParams(params);
        Vector modeled = _evaluator->computeModelImage();
        Vector residual = (_measured-modeled);
        Matrix const & lpd = _evaluator->computeLinearParameterDerivative();
        Matrix const & npd = _evaluator->computeNonlinearParameterDerivative();

        Eigen::VectorXd priorDelta(params.size());
        std::copy(params.begin(), params.end(), priorDelta.data());
        priorDelta -= _priorMean;

        std::vector<double> gradient(params.size());
        multifit::VectorMap gradMap(&gradient[0], params.size());
        gradMap.start(lpd.cols()) = (-lpd).transpose()*residual;
        gradMap.end(npd.cols()) = (-npd).transpose()*residual;
        gradMap += _priorFisherDiag.cwise() * priorDelta;
        return gradient; 
    }
private:
    bool _dirty;
    bool _doSnapshots;
    int _nUpdates;
    std::string _snapshotFormat;
    multifit::ModelEvaluator::Ptr _evaluator;
    Vector _measured;
    Vector _priorMean;
    Vector _priorFisherDiag;

    void checkParams(std::vector<double> const & params) {
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

    void setParams(std::vector<double> const & params) {
        checkParams(params);
        if(_dirty) {
            int offset = 0;
            if(_evaluator->getLinearParameterSize() > 0) {
                _evaluator->setLinearParameters(&params.front());
                offset += _evaluator->getLinearParameterSize();
            }
            if(_evaluator->getNonlinearParameterSize() > 0 ) {
                _evaluator->setNonlinearParameters(&params.front() + offset);
            }
        }
        if (_doSnapshots) {
            _evaluator->getProjectionList().front()->writeSnapshot(
                boost::str(boost::format(_snapshotFormat) % _nUpdates),
                _evaluator->getDataVector()
            );
        ++_nUpdates;
        }
    }

};


class Function : public ROOT::Minuit2::FCNBase {
public:
    Function(ChisqFunction const & chisqFunction) :
        _chisqFunction(chisqFunction)
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
        ChisqFunction const & chisqFunction,
        bool checkGradient=false
    ) :
        _chisqFunction(chisqFunction),
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


multifit::MinuitFitterResult doMinuitFit(
    multifit::ModelEvaluator::Ptr evaluator, 
    std::vector<double> const & initialErrors,
    std::vector<double> const & priorMean,
    std::vector<double> const & priorFisherDiag,
    std::vector<double> const & lower,
    std::vector<double> const & upper,
    lsst::pex::policy::Policy::Ptr policy
) {
    int nLinear = evaluator->getLinearParameterSize();
    int nNonlinear = evaluator->getNonlinearParameterSize();
    int nParams = nLinear + nNonlinear;

    if(nParams != int(initialErrors.size())) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Number of model parameters not equal to length of error vector"
        );
    }
    if(nParams != int(priorMean.size()) || nParams != int(priorFisherDiag.size())) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Number of model parameters not equal to length of prior mean and/or Fisher diagonal vector"
        );
    }

    if(nParams != int(lower.size()) || nParams != int(upper.size())) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Number of model parameters not equal to length of limits vectors"
        );
    }

    ChisqFunction::Vector priorMeanEig(priorMean.size());
    std::copy(priorMean.begin(), priorMean.end(), priorMeanEig.data());

    ChisqFunction::Vector priorFisherDiagEig(priorFisherDiag.size());
    std::copy(priorFisherDiag.begin(), priorFisherDiag.end(), priorFisherDiagEig.data());

    ChisqFunction chisqFunction(
        evaluator, priorMeanEig, priorFisherDiagEig, 
        policy->getBool("doSnapshots"),
        policy->getString("snapshotFormat")
    );

    std::vector<double> initialParams(nParams);
    multifit::VectorMap paramMap(&initialParams[0], nParams);
    if(nLinear > 0){
        paramMap.start(nLinear) = evaluator->getLinearParameters();
    }
    if(nNonlinear > 0){
        paramMap.end(nNonlinear) = evaluator->getNonlinearParameters();
    }
    ROOT::Minuit2::MnUserParameters userParams(initialParams, initialErrors);
    for(int i =0; i < nParams; ++i) {
        double const & l = lower[i];
        double const & u = upper[i];
        if(l != std::numeric_limits<double>::infinity()){
            if(u != std::numeric_limits<double>::infinity()) {
                userParams.SetLimits(i, l, u);
            }
            else {
                userParams.SetLowerLimit(i, l);
            }
        }
        else if(u != std::numeric_limits<double>::infinity()) {
            userParams.SetUpperLimit(i, u);
        }
    }

    boost::scoped_ptr<ROOT::Minuit2::FCNBase> function;
    boost::scoped_ptr<ROOT::Minuit2::MnApplication> minimizer;
    std::string algorithm = policy->getString("algorithm"); 
    int strategy = policy->getInt("strategy");
    if(algorithm == "MIGRAD") {
        if(policy->getBool("doAnalyticGradient")) {
            bool checkGradient = policy->getBool("checkGradient");
            function.reset(new GradientFunction(chisqFunction, checkGradient));
            minimizer.reset(
                new ROOT::Minuit2::MnMigrad(
                    static_cast<ROOT::Minuit2::FCNGradientBase&>(*function), userParams, strategy
                )
            );
        } else {
            function.reset(new Function(chisqFunction));
            minimizer.reset(new ROOT::Minuit2::MnMigrad(*function, userParams, strategy));
        }
    } else if (algorithm == "SIMPLEX") {
        function.reset(new Function(chisqFunction));
        minimizer.reset(new ROOT::Minuit2::MnSimplex(*function, userParams, strategy));
    }

    ROOT::Minuit2::FunctionMinimum min = (*minimizer)(
        policy->getInt("iterationMax"), 
        policy->getDouble("tolerance")
    );
    return multifit::MinuitFitterResult(min);
}

}//end anonymous namespace


/******************************************************************************
 * MinuitFitter
 *****************************************************************************/
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

multifit::MinuitFitterResult multifit::MinuitFitter::apply(
    multifit::ModelEvaluator::Ptr evaluator, 
    std::vector<double> const & initialErrors
) const {
    int nParams = evaluator->getLinearParameterSize() + evaluator->getNonlinearParameterSize();
    std::vector<double> priorMean(nParams, 0.0);
    std::vector<double> priorFisherDiag(nParams, 0.0);
    return apply(evaluator, initialErrors, priorMean, priorFisherDiag);
}

multifit::MinuitFitterResult multifit::MinuitFitter::apply(
    multifit::ModelEvaluator::Ptr evaluator, 
    std::vector<double> const & initialErrors,
    std::vector<double> const & priorMean,
    std::vector<double> const & priorFisherDiag
) const {
    int nParams = evaluator->getLinearParameterSize() + evaluator->getNonlinearParameterSize();
    std::vector<double> lower(nParams, std::numeric_limits<double>::infinity());
    std::vector<double> upper(nParams, std::numeric_limits<double>::infinity());
    return apply(evaluator, initialErrors, priorMean, priorFisherDiag, lower, upper);
}

multifit::MinuitFitterResult multifit::MinuitFitter::apply(
    multifit::ModelEvaluator::Ptr evaluator, 
    std::vector<double> const & initialErrors,
    std::vector<double> const & priorMean,
    std::vector<double> const & priorFisherDiag,
    std::vector<double> const & lowerLimits,
    std::vector<double> const & upperLimits
) const {
    return ::doMinuitFit(
        evaluator, initialErrors, 
        priorMean, priorFisherDiag, 
        lowerLimits, upperLimits, 
        _policy
    ); 
}
