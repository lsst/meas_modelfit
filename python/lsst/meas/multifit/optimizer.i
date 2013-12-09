%{
#include "lsst/meas/multifit/optimizer.h"
%}

%returnCopy(Optimizer::getControl)
%returnCopy(Optimizer::getIterations)
%shared_ptr(lsst::meas::multifit::OptimizerObjective)

%include "lsst/meas/multifit/optimizer.h"

%usePointerEquality(lsst::meas::multifit::OptimizerObjective)

// ignore std::vector methods that require a default ctor, as we don't have one for this class
%ignore std::vector<lsst::meas::multifit::OptimizerIterationData>::vector(size_type);
%ignore std::vector<lsst::meas::multifit::OptimizerIterationData>::resize(size_type);
%template(OptimizerIterationDataVector) std::vector<lsst::meas::multifit::OptimizerIterationData>;

%extend lsst::meas::multifit::Optimizer {
    %pythoncode %{
        IterationData = OptimizerIterationData
        Objective = OptimizerObjective
        Control = OptimizerControl
        HistoryRecorder = OptimizerHistoryRecorder

        def getConfig(self):
            config = self.ConfigClass()
            config.readControl(self.getControl())
            return config
    %}
}

%pythoncode %{
    import lsst.pex.config
    @lsst.pex.config.wrap(OptimizerControl)
    class OptimizerConfig(lsst.pex.config.Config):
        pass
    Optimizer.ConfigClass = OptimizerConfig
    Optimizer.IterationDataVector = OptimizerIterationDataVector
%}
