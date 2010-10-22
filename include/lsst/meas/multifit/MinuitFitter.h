#ifndef LSST_MEAS_MULTIFIT_MINUITFITTER_H
#define LSST_MEAS_MULTIFIT_MINUITFITTER_H

#include <Minuit2/FunctionMinimum.h>
#include <float.h>
#include "lsst/pex/policy/DefaultPolicyFile.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/meas/multifit/ModelEvaluator.h"

namespace lsst {
namespace meas {
namespace multifit {

struct MinuitFitterResult {
public:
    enum Flags {
        CONVERGED = 1,
        MAX_ITERATION_REACHED = 2
    };
    
#ifndef SWIG
    MinuitFitterResult(ROOT::Minuit2::FunctionMinimum const & min) 
      : flags(0),  chisq(min.Fval()), nIterations(min.NFcn()) {
        if(min.IsValid()) 
            flags |= CONVERGED;
        if(min.HasReachedCallLimit())
            flags |= MAX_ITERATION_REACHED;
        parameters = min.UserParameters().Params();
        int nParam=parameters.size();
        covariance = Eigen::MatrixXd(nParam, nParam);
        for (int i=0; i < nParam; ++i) {
            for (int j=0; j< nParam; ++j) {
                covariance(i,j) = min.UserCovariance()(i,j);
            }
        }
    };
#endif
    
    MinuitFitterResult() : flags(0), chisq(DBL_MAX), nIterations(0), model() {}

    int flags;
    double chisq;
    int nIterations;
    std::vector<double> parameters;
    Model::ConstPtr model;
    Eigen::MatrixXd covariance;
};

class MinuitAnalyticFitter {
public:

    typedef MinuitFitterResult Result;

    MinuitAnalyticFitter(
        lsst::pex::policy::Policy::Ptr const & policy = lsst::pex::policy::Policy::Ptr()
    );
    
    Result apply(multifit::ModelEvaluator::Ptr evaluator, std::vector<double> const & initialErrors) const;

    Result apply(
        multifit::ModelEvaluator::Ptr evaluator,
        std::vector<double> const & initialErrors,
        std::vector<double> const & priorMean,
        std::vector<double> const & priorFisherDiag
    ) const;

    /**
     * Retrieve the default policy for configuring a SingleLinearParameterFitter
     *
     * The defaults are provided in the dictionary file
     */
    static PTR(lsst::pex::policy::PolicySource) getDefaultPolicySource() {
        static const PTR(lsst::pex::policy::PolicySource) source(
            new lsst::pex::policy::DefaultPolicyFile(
                "meas_multifit", 
                "MinuitFitterDict.paf", 
                "policy"
            )
        );
        return source;
    }
private:
    lsst::pex::policy::Policy::Ptr _policy;
};

class MinuitNumericFitter {
public:

    typedef MinuitFitterResult Result;

    MinuitNumericFitter(lsst::pex::policy::Policy::Ptr const & policy = lsst::pex::policy::Policy::Ptr());
    
    Result apply(multifit::ModelEvaluator::Ptr evaluator, std::vector<double> const & initialErrors) const;

    Result apply(
        multifit::ModelEvaluator::Ptr evaluator,
        std::vector<double> const & initialErrors,
        std::vector<double> const & priorMean,
        std::vector<double> const & priorFisherDiag
    ) const;

    /**
     * Retrieve the default policy for configuring a SingleLinearParameterFitter
     *
     * The defaults are provided in the dictionary file
     */
    static PTR(lsst::pex::policy::PolicySource) getDefaultPolicySource() {
        static const PTR(lsst::pex::policy::PolicySource) source(
            new lsst::pex::policy::DefaultPolicyFile(
                "meas_multifit", 
                "MinuitFitterDict.paf", 
                "policy"
            )
        );
        return source;
    }
private:
    lsst::pex::policy::Policy::Ptr _policy;
};

}}}

#endif
