#include "lsst/meas/multifit/LevMarFitter.h"
#include "levmar/levmar.h"
#include <Eigen/Core>
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/pex/logging/Trace.h"
#include <limits>
#include <iostream>
#include <boost/format.hpp>
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

    explicit LevMarFunction(
        multifit::ModelEvaluator::Ptr const & evaluator,
        std::string const & snapshotFormat, bool doSnapshots=false
    ) 
        : _dirty(true), _doSnapshots(doSnapshots), 
          _nUpdates(0), _snapshotFormat(snapshotFormat), _evaluator(evaluator)
    {}

    void computeModel(VectorMap const & params, VectorMap & model) {
        setParams(params);
        model = _evaluator->computeModelImage();
    }

    void computeJacobian(VectorMap const & params, MatrixMap & jacobian) {
        setParams(params);
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
        self->computeModel(params, model);
    }

    static void jacf(double * p, double * j, int m, int n, void * data) {
        VectorMap params(p, m);
        MatrixMap jacobian(j, n, m);
        LevMarFunction * self = reinterpret_cast<LevMarFunction*>(data);
        self->computeJacobian(params, jacobian);
    }

private:
    bool _dirty;
    bool _doSnapshots;
    int _nUpdates;
    std::string _snapshotFormat;
    multifit::ModelEvaluator::Ptr _evaluator;

    void checkParams(VectorMap const & params) {
        if(_dirty)
            return;
        
        double const * iNew = params.data();
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

    void setParams(VectorMap const & params) {
        checkParams(params);
        if(_dirty) {
            std::cerr << "stepping to params: " << params << std::endl;
            if(_evaluator->getLinearParameterSize() > 0) {
                _evaluator->setLinearParameters(
                    params.start(_evaluator->getLinearParameterSize())
                );
            }
            if(_evaluator->getNonlinearParameterSize() > 0 ) {
                _evaluator->setNonlinearParameters(
                    params.end(_evaluator->getNonlinearParameterSize())
                );
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

}//end anonymous namespace

multifit::LevMarFitterResult::LevMarFitterResult(
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
    if(nLinearParams > 0) {
        parameters.start(nLinearParams) = evaluator->getLinearParameters();
    }
    if(nNonlinearParams > 0) {
        parameters.end(nNonlinearParams) = evaluator->getNonlinearParameters();
    }
    Eigen::VectorXd result(evaluator->getNPixels());
    LevMarFunction adata(evaluator, _policy->getString("snapshotFormat"), _policy->getBool("doSnapshots"));
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
    if(nLinearParams > 0)
        parameters.start(nLinearParams) = evaluator->getLinearParameters();
    if(nNonlinearParams > 0)
        parameters.end(nNonlinearParams) = evaluator->getNonlinearParameters();

    LevMarFunction::MatrixRM covar(parameters.size(), parameters.size());
    LevMarFunction::Vector image = evaluator->getWeightedData();
    LevMarFunction adata(evaluator, _policy->getString("snapshotFormat"), _policy->getBool("doSnapshots"));
    double info[LM_INFO_SZ];
    int itmax = _policy->getInt("iterationMax");
    double opts[5] = {
        _policy->getDouble("tau"),
        _policy->getDouble("gradientEpsilon"),
        _policy->getDouble("parameterEpsilon"),
        _policy->getDouble("residualEpsilon")
    };
    if (_policy->getBool("doAnalyticJacobian")) {
        (void)dlevmar_der(
            &LevMarFunction::func, &LevMarFunction::jacf, parameters.data(), image.data(), 
            parameters.size(), image.size(), itmax, opts, info, 0, covar.data(), &adata
        );
    } else {
        opts[4] = _policy->getDouble("delta");
        (void)dlevmar_dif(
            &LevMarFunction::func, parameters.data(), image.data(), 
            parameters.size(), image.size(), itmax, opts, info, 0, covar.data(), &adata
        );
    }
    return Result(info, evaluator->getModel(), parameters, covar);
}
