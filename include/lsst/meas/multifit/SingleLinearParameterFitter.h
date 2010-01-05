#ifndef LSST_MEAS_MULTIFIT_SINGLE_LINEAR_PARAMETER_FITTER_H
#define LSST_MEAS_MULTIFIT_SINGLE_LINEAR_PARAMETER_FITTER_H

#include "lsst/pex/policy/PolicyConfigured.h"
#include "lsst/daf/base/PropertySet.h"

#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/ModelEvaluator.h"
#include "lsst/meas/multifit/matrices.h"

namespace lsst {
namespace meas {
namespace multifit {
class SingleLinearParameterFitter : public lsst::pex::policy::PolicyConfigured {
public:
    enum TerminationType {ITERATION=1, CHISQ=2, STEP=4};

    struct Result {
        typedef boost::shared_ptr<Result> Ptr;

        Result() : convergenceFlags(0), chisq(0), model(), sdqaMetrics() {}
        int convergenceFlags;
        double chisq, dChisq;
        Model::ConstPtr model;
        lsst::daf::base::PropertySet sdqaMetrics;
    };

    SingleLinearParameterFitter(lsst::pex::policy::Policy::Ptr const & policy);

    Result::Ptr apply(ModelEvaluator &) const;

private:
    bool terminateIteration(
        double const & chisq, int const & nIterations, Eigen::VectorXd const & step
    ) const {
        return ((terminationType & ITERATION) && nIterations >= iterationMax) ||
            ((terminationType & CHISQ) && chisq <= chisqThreshold) ||
            ((terminationType & STEP) && step.norm() <= stepThreshold); 
     }

    int terminationType;
    int iterationMax;
    double chisqThreshold;
    double stepThreshold;
    
};

}}} //end namespace lsst::meas::multifit

#endif //LSST_MEAS_MULTIFIT_SINGLE_LINEAR_PARAMETER_FITTER_H
