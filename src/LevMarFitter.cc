#include "lsst/meas/multifit/LevMarFitter.h"
#include "levmar/levmar.h"
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

class LevMarFunction {
public:
    typedef Eigen::Matrix<multifit::Pixel, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Map<Vector> VectorMap;
    typedef Eigen::Matrix<multifit::Pixel, Eigen::Dynamic, Eigen::Dynamic,Eigen::RowMajor> MatrixRM;
    typedef Eigen::Matrix<multifit::Pixel, Eigen::Dynamic, Eigen::Dynamic> MatrixCM;
    typedef Eigen::Map<MatrixRM> MatrixMap;

    explicit LevMarFunction(multifit::ModelEvaluator::Ptr const & evaluator) 
      : _dirty(true), 
        _evaluator(evaluator)
    {}

    void computeModel(VectorMap const & params, VectorMap & model) {
        checkParams(params);
        if (_dirty) {
            _evaluator->setLinearParameters(&params[0]);
            _evaluator->setNonlinearParameters(&params[_evaluator->getLinearParameterSize()]);
        }
        model = _evaluator->computeModelImage();
    }

    void computeJacobian(VectorMap const & params, MatrixMap & jacobian) {
        checkParams(params);
        if(_dirty) {
            _evaluator->setLinearParameters(&params[0]);
            _evaluator->setNonlinearParameters(&params[_evaluator->getLinearParameterSize()]);            
        }
        MatrixCM const & lpd = _evaluator->computeLinearParameterDerivative();
        MatrixCM const & npd = _evaluator->computeNonlinearParameterDerivative();
        assert(jacobian.cols() == lpd.cols() + npd.cols());
        jacobian.block(0, 0, jacobian.rows(), lpd.cols()) = lpd;
        jacobian.block(0, lpd.cols(), jacobian.rows(), npd.cols()) = npd;
    }

    static void func(double * p, double * hx, int m, int n, void * data) {
        VectorMap params(p, m);
        VectorMap model(hx, n);
        LevMarFunction * self = reinterpret_cast<LevMarFunction*>(data);
        self->computeModel(params, residuals);
    }

    static void jacf(double * p, double * j, int m, int n, void * data) {
        VectorMap params(p, m);
        MatrixMap jacobian(j, n, m);
        LevMarFunction * self = reinterpret_cast<LevMarFunction*>(data);
        self->computeJacobian(params, jacobian);
    }

private:
    bool _dirty;
    multifit::ModelEvaluator::Ptr _evaluator;

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

}//end anonymous namespace

LevMarFitterResult::LevMarFitterResult(
    double * levmarInfo,
    Model::ConstPtr const & model_,
    Eigen::VectorXd const & parameters_,
    Eigen::MatrixXd const & covariance_
) :
    termination(TerminationEnum(int(levmarInfo[6]))),
    chisqInitial(levmarInfo[0]),
    chisqFinal(levmarInfo[1]),
    maxGradient(levmarInfo[2]),
    lastStepNorm(levmarInfo[3]),
    nIterations(int(levmarInfo[5])),
    nFunctionEvaluations(levmarInfo[7]),
    nJacobianEvaluations(levmarInfo[8]),
    nMatrixFactorizations(levmarInfo[9]),
    model(model_),
    parameters(parameters_),
    covariance(covariance_)
{}

/******************************************************************************
 * LevMarFitter
 *****************************************************************************/
multifit::LevMarFitter::LevMarFitter(
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

Eigen::VectorXd multifit::LevMarFitter::checkDerivatives(
    multifit::ModelEvaluator::Ptr const & evaluator
) const {
    int nLinearParams = evaluator->getLinearParameterSize();
    int nNonlinearParams = evaluator->getNonlinearParameterSize();
    LevMarFunction::Vector parameters(nLinearParams + nNonlinearParams);
    parameters.segment(0, nLinearParams) = evaluator->getLinearParameters();
    parameters.segment(nLinearParams, nNonLinearParams) = evaluator->getNonlinearParameters();
    Eigen::VectorXd result(evaluator->getNPixels());
    LevMarFunction adata(evaluator);
    dlevmar_chkjac(
        &LevMarFunction::func, &LevMarFunction::jacf, parameters.data(), 
        parameters.size(), result.size(), &adata, result.data()
    );
    return result;
}


multifit::LevMarFitterResult multifit::LevMarFitter::apply(
    multifit::ModelEvaluator::Ptr const & evaluator
) const {
    int nLinearParams = evaluator->getLinearParameterSize();
    int nNonlinearParams = evaluator->getNonlinearParameterSize();
    LevMarFunction::Vector parameters(nLinearParams + nNonlinearParams);
    parameters.segment(0, nLinearParams) = evaluator->getLinearParameters();
    parameters.segment(nLinearParams, nNonLinearParams) = evaluator->getNonlinearParameters();
    LevMarFunction::MatrixRM covariance(parameters.size(), parameters.size());
    LevMarFunction::Vector image = evaluator->getWeightedData();
    LevMarFunction adata(evaluator);
    double info[LM_INFO_SZ];
    int itmax = policy->getInt("iterationMax");
    int nIterations = 0;
    double opts[5] = {
        policy->getDouble("tau"),
        policy->getDouble("gradientEpsilon"),
        policy->getDouble("parameterEpsilon"),
        policy->getDouble("residualEpsilon")
    };
    if (policy->getBool("doAnalyticJacobian")) {
        nIterations = dlevmar_der(
            &LevMarFunction::func, &LevMarFunciton::jacf, parameters.data(), image.data(), 
            parameter.size(), image.size(), itmax, opts, info, 0, covar.data(), &adata
        );
    } else {
        opts[4] = policy->getDouble("delta");
        nIterations = dlevmar_dif(
            &LevMarFunction::func, parameters.data(), image.data(), 
            parameter.size(), image.size(), itmax, opts, info, 0, covar.data(), &adata
        );
    }
    return Result(info, nIterations, parameters, covariance);
}
