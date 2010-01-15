#ifndef LSST_MEAS_MULTIFIT_SINGLE_LINEAR_PARAMETER_FITTER_H
#define LSST_MEAS_MULTIFIT_SINGLE_LINEAR_PARAMETER_FITTER_H

#include <float.h>

#include "lsst/pex/policy/PolicyConfigured.h"
#include "lsst/daf/base/PropertySet.h"

#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/ModelEvaluator.h"
#include "lsst/meas/multifit/matrices.h"

namespace lsst {
namespace meas {
namespace multifit {



struct SimpleFitResult {
    enum BasicConvergenceFlags {
        CONVERGED = 1,
        MAX_ITERATION_REACHED = 2,
        DCHISQ_THRESHOLD_REACHED = 4,
        STEP_THRESHOLD_REACHED = 8
    };

    typedef boost::shared_ptr<SimpleFitResult> Ptr;

    SimpleFitResult() 
        : convergenceFlags(0), 
          chisq(DBL_MAX), dChisq(DBL_MAX),            
          sdqaMetrics(new lsst::daf::base::PropertySet())
    {}

    int convergenceFlags;
    double chisq, dChisq;
    Model::ConstPtr model;
    lsst::daf::base::PropertySet::Ptr sdqaMetrics;
};

class SingleLinearParameterFitter : public lsst::pex::policy::PolicyConfigured {
public:
    enum TerminationType {ITERATION=1, DCHISQ=2, STEP=4};

    typedef SimpleFitResult Result;

    SingleLinearParameterFitter(lsst::pex::policy::Policy::Ptr const & policy);

    Result::Ptr apply(ModelEvaluator &) const;

private:
    int _terminationType;
    int _iterationMax;
    double _dChisqThreshold;
    double _stepThreshold;
    
};

}}} //end namespace lsst::meas::multifit

#endif //LSST_MEAS_MULTIFIT_SINGLE_LINEAR_PARAMETER_FITTER_H
