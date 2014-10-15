%{
#include "lsst/meas/modelfit/optimizer.h"
%}

%returnCopy(Optimizer::getControl)
%returnCopy(Optimizer::getIterations)
%shared_ptr(lsst::meas::modelfit::OptimizerObjective)
%shared_ptr(lsst::meas::modelfit::OptimizerInterpreter)

%include "lsst/meas/modelfit/optimizer.h"

%castShared(lsst::meas::modelfit::OptimizerInterpreter, lsst::meas::modelfit::Interpreter)
%usePointerEquality(lsst::meas::modelfit::OptimizerObjective)

// ignore std::vector methods that require a default ctor, as we don't have one for this class
%ignore std::vector<lsst::meas::modelfit::OptimizerIterationData>::vector(size_type);
%ignore std::vector<lsst::meas::modelfit::OptimizerIterationData>::resize(size_type);
%template(OptimizerIterationDataVector) std::vector<lsst::meas::modelfit::OptimizerIterationData>;

%extend lsst::meas::modelfit::Optimizer {
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
