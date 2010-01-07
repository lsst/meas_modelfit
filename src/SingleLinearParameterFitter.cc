#include <float.h>

#include <Eigen/Core>
#include <Eigen/QR>
#include <Eigen/Cholesky>

#include "lsst/meas/multifit/SingleLinearParameterFitter.h"
#include "lsst/pex/exceptions/Runtime.h"

namespace multifit = lsst::meas::multifit;

multifit::SingleLinearParameterFitter::SingleLinearParameterFitter(
    lsst::pex::policy::Policy::Ptr const & policy
) : lsst::pex::policy::PolicyConfigured(policy) {
    if(!policy->exists("terminationType")) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Invalid configuration policy - missing value \"terminationType\""
        );
    } 
    
    terminationType = 0;
    std::vector<std::string> terminationVector(policy->getStringArray("terminationType"));
    for(std::vector<std::string>::iterator i(terminationVector.begin()), end(terminationVector.end());
        i != end; ++i
    ) {
        if((*i) == "chisq")
            terminationType |= CHISQ;
        else if((*i) == "iteration")
            terminationType |= ITERATION;
        else if((*i) == "step")
            terminationType |= STEP; 
    }
    
    if(terminationType & ITERATION) {
        if(!policy->exists("iterationMax")) {
            throw LSST_EXCEPT(
                lsst::pex::exceptions::InvalidParameterException,
                "Invalid configuration policy - missing value \"iterationMax\""
            );
        } else {
            iterationMax = policy->getInt("iterationMax");   
        }
    }
    if(terminationType & CHISQ) {
        if(!policy->exists("chisqThreshold")) {
            throw LSST_EXCEPT(
                lsst::pex::exceptions::InvalidParameterException,
                "Invalid configuration policy - missing value \"chisqThreshold\""
            );
        } else {
            chisqThreshold = policy->getDouble("chisqThreshold");   
        }
    }
    if(terminationType & STEP) {
        if(!policy->exists("stepThreshold")) {
            throw LSST_EXCEPT(
                lsst::pex::exceptions::InvalidParameterException,
                "Invalid configuration policy - missing value \"stepThreshold\""
            );
        } else {
            stepThreshold = policy->getDouble("stepThreshold");   
        }
    }

    lsst::pex::policy::PolicyConfigured::configured();
}    

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
    int nIterations = 0;

    VectorMap image(
        evaluator.getImageVector().getData(),
        evaluator.getNPixels()
    );
    VectorMap variance(
        evaluator.getVarianceVector().getData(),
        evaluator.getNPixels()
    );
    if(!(image.size() >0 && variance.size()> 0))
        throw "foobared, something went to hell";

    Eigen::VectorXd data = (image.cwise() / (variance.cwise().sqrt()));

    Eigen::VectorXd newLinearDNonlinear;
    Eigen::MatrixXd jacobian;

    double chisq, dChisq=DBL_MAX; 
    Eigen::VectorXd step;
    do {        
        VectorMap dLinear = VectorMap(
            evaluator.computeLinearParameterDerivative().getData(),
            evaluator.getNPixels()
        );
        MatrixMap dNonlinear = MatrixMap(
            evaluator.computeNonlinearParameterDerivative().getData(),
            evaluator.getNPixels(),
            evaluator.getNonlinearParameterSize()
        );


        Eigen::QR<Eigen::VectorXd> qrLinear(dLinear);
        
        //qr.matrixR() is a matrix with dimensions (1, 1), represent it as a double
        double R = qrLinear.matrixR()(0,0);

        //the new linear parameter as a funtion of the nonlinear parameters
        double newLinear = (qrLinear.matrixQ().transpose()*data)(0,0)/R;
   
        //compute the residual as a function of the nonlinear parameters
        Eigen::VectorXd residual = data - dLinear*newLinear;

        //compute the chisq
        dChisq = chisq;
        chisq = (residual.dot(residual))/2;
        dChisq -= chisq;
        

        //compute derivative of new linear parameter w.r.t nonlinear parameters
        //this is a matrix with dimensions (nNonlinear, 1)
        Eigen::VectorXd dNewLinear = dNonlinear.transpose()*(residual - dLinear*newLinear);
        dNewLinear /= (R*R*newLinear);
    
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
    } while(!terminateIteration(dChisq, nIterations, step));

    result->chisq = chisq;
    result->dChisq = dChisq;
    result->sdqaMetrics->set("nIterations", nIterations);

    return result;
}

