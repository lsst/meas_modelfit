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

    terminationType = policy->getInt("terminationType");
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
    Result::Ptr result = boost::make_shared<Result>();
    result->model = evaluator.getModel();    
    int nIterations = 0;
    double chisq;

    VectorMap image = getVectorView(evaluator.getImageVector());
    VectorMap variance = getVectorView(evaluator.getVarianceVector());
    Eigen::VectorXd data = image.cwise() / variance.cwise().sqrt();

    Eigen::VectorXd linear;
    Eigen::VectorXd nonlinear;

    Eigen::VectorXd residual;
    double newLinear;
    Eigen::VectorXd newLinearDNonlinear;
    Eigen::MatrixXd jacobian;
    
    Eigen::VectorXd step;
    do {        
        MatrixMap dLinear = MatrixMap(
            evaluator.computeLinearParameterDerivative().getData(),
            evaluator.getNPixels(),
            evaluator.getLinearParameterSize()
        );
        MatrixMap dNonlinear = MatrixMap(
            evaluator.computeNonlinearParameterDerivative().getData(),
            evaluator.getNPixels(),
            evaluator.getNonlinearParameterSize()
        );


        Eigen::QR<Eigen::MatrixXd> qrLinear(dLinear);
        
        //qr.matrixR() is a matrix with dimensions (1, 1), represent it as a double
        double R = qrLinear.matrixR()(0,0);

        //the new linear parameter as a funtion of the nonlinear parameters
        newLinear = (qrLinear.matrixQ().transpose()*data)(0,0)/R;
   
        //compute the residual as a function of the nonlinear parameters
        residual = data - dLinear*newLinear;

        //compute the chisq
        chisq = (residual.dot(residual))/2;

        //compute derivative of new linear parameter w.r.t nonlinear parameters
        //this is a matrix with dimensions (nNonlinear, 1)
        newLinearDNonlinear = dNonlinear.transpose()*(residual - dLinear*newLinear);
        newLinearDNonlinear /= (R*R*newLinear);
    
        //compute the jacobian of partial derivatives of the model w.r.t
        //nonlinear parameters
        //this is a matrix with dimensions (pixels, nNonlinear)
        jacobian = -dNonlinear - dLinear*newLinearDNonlinear.transpose();

        //compute the step to take on nonlinear parameters:
        //this is a matrix with dimensions (nNonlinear, 1)
        step = jacobian.transpose()*residual;

        (jacobian.transpose()*jacobian).llt().solveInPlace(step);

        linear[0] = newLinear;
        nonlinear = evaluator.getNonlinearParameters() + step;
        evaluator.setLinearParameters(linear.data());
        evaluator.setNonlinearParameters(nonlinear.data());
    } while(!terminateIteration(chisq, nIterations, step));

    result->chisq = chisq;
    result->sdqaMetrics.set("nIterations", nIterations);

    return result;
}

