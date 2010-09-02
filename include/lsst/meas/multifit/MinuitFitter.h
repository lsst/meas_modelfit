#ifndef LSST_MEAS_MULTIFIT_MINUITFITTER_H
#define LSST_MEAS_MULTIFIT_MINUITFITTER_H

#include <Minuit2/FunctionMinimum.h>
#include "lsst/pex/policy/DefaultPolicyFile.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/meas/multifit/ModelEvaluator.h"

namespace lsst {
namespace meas {
namespace multifit {

class MinuitFitter {
public:
    enum ConvergenceFlags {
        CONVERGED = 1,
        MAX_ITERATION_REACHED = 2
    };

    struct Result {
    public:
        Result(ROOT::Minuit2::FunctionMinimum const & min) : _min(min), _flags(0){
            if(_min.IsValid()) 
                _flags |= CONVERGED;
            if(_min.HasReachedCallLimit())
                _flags |= MAX_ITERATION_REACHED;
        };
        int getConvergenceFlags() const {return _flags;}
        double getChisq() const {return _min.Fval();}
        double getNIterations() const {return _min.NFcn();}
    private:
        ROOT::Minuit2::FunctionMinimum _min;
        int _flags;
    };

    MinuitFitter(lsst::pex::policy::Policy::Ptr const & policy);
    
    Result apply(multifit::ModelEvaluator::Ptr evaluator, std::vector<double> initialErrors) const;
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
