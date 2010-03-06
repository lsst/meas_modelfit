// -*- lsst-c++ -*-
/**
 * Support for fitting models with a single linear parameter
 *
 * Contains declarions for SimpleFitResult and SingleLinearParameterFitter
 */
#ifndef LSST_MEAS_MULTIFIT_SINGLE_LINEAR_PARAMETER_FITTER_H
#define LSST_MEAS_MULTIFIT_SINGLE_LINEAR_PARAMETER_FITTER_H

#include <float.h>


#include "lsst/pex/policy/DefaultPolicyFile.h"
#include "lsst/daf/base/PropertySet.h"

#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/ModelEvaluator.h"

namespace lsst {
namespace meas {
namespace multifit {


/**
 * Bare-bones model fitting result
 *
 * Records the state upon termination of the model-fitting loop. This includes
 * information about if and how convergence was determined, as well as the
 * chisq value, chisq delta of the last iteration, and a property which fitters
 * may use to hold aditional metrics
 */
struct SimpleFitResult {
    enum ConvergenceFlags {
        CONVERGED = 1,
        MAX_ITERATION_REACHED = 2,
        DCHISQ_THRESHOLD_REACHED = 4,
        STEP_THRESHOLD_REACHED = 8
    };

    typedef boost::shared_ptr<SimpleFitResult> Ptr;

    /**
     * Default construct a SimpleFitResult
     */
    SimpleFitResult() 
        : convergenceFlags(0), 
          chisq(DBL_MAX), dChisq(DBL_MAX),            
          sdqaMetrics(new lsst::daf::base::PropertySet())
    {}
    /// bit flags indicating if and how the fitter converged
    int convergenceFlags;     
    /// chisq value of the last iteration
    double chisq;
    /// delta chisq of between the last two iterations
    double dChisq;
    /// pointer to the model that was fit
    Model::ConstPtr model;
    /// additional properties set by the fitter
    lsst::daf::base::PropertySet::Ptr sdqaMetrics;
};

/**
 * A model fitter optimized for a single linear parameter
 *
 * This fitter is a slight modification on the algorithm of seperable nonlinear
 * least-squares, which takes into account that there is only one linear
 * parameter.
 *
 * The fitter can be configured via policy upon construction. A policy
 * dictionary (with default values) can be found at:
 * lsst/meas/multifit/policy/SingleLinearParameterFitterDict.paf
 *
 * Three termination conditions are supported:
 * @li \c iteration once a certain number of iterations are performed, the fitter
 *      stops. This is useful for debugging, and for setting an upper-limit on
 *      runtime, if worried that a model will not converge.
 * @li \c dChisq Consider fitter to converge when the difference in chisq
 *      between two iterations drops below threshold.
 * @li \c step Consider fitter to converge when the norm of the step in nonlinear 
 *      parameters drops below threshold.
 */
class SingleLinearParameterFitter  {
public:
    typedef SimpleFitResult Result;
    typedef lsst::pex::policy::Policy::Ptr PolicyPtr;
    typedef boost::shared_ptr<lsst::pex::policy::PolicySource> PolicySourcePtr;

    enum TerminationType {ITERATION=1, DCHISQ=2, STEP=4};

    SingleLinearParameterFitter(PolicyPtr const & policy=PolicyPtr());

    Result::Ptr apply(ModelEvaluator &) const;

    /**
     * Retrieve the default policy for configuring a SingleLinearParameterFitter
     *
     * The defaults are provided in the dictionary file
     */
    static PolicySourcePtr getDefaultPolicySource() {
        static const PolicySourcePtr source(
            new lsst::pex::policy::DefaultPolicyFile(
                "meas_multifit", 
                "SingleLinearParameterFitterDict.paf", 
                "policy"
            )
        );
        return source;
    }
private:
    int _terminationType;
    int _iterationMax;
    double _dChisqThreshold;
    double _stepThreshold;
    PolicyPtr _policy; 
};

}}} //end namespace lsst::meas::multifit

#endif //LSST_MEAS_MULTIFIT_SINGLE_LINEAR_PARAMETER_FITTER_H
