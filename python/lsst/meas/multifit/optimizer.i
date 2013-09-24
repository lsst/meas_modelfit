%{
#include "lsst/meas/multifit/optimizer.h"
%}

%declareNumPyConverters(ndarray::Array<double,1,1>)
%declareNumPyConverters(ndarray::Array<double const,1,1>)
%declareNumPyConverters(ndarray::Array<double,2,1>)
%declareNumPyConverters(ndarray::Array<double const,2,1>)
%declareNumPyConverters(ndarray::Array<double,2,2>)
%declareNumPyConverters(ndarray::Array<double const,2,2>)

%returnCopy(PosteriorOptimizer::getControl)
%returnCopy(PosteriorOptimizer::getIterations)
%shared_ptr(lsst::meas::multifit::PosteriorOptimizerObjective)

%include "lsst/meas/multifit/optimizer.h"

%usePointerEquality(lsst::meas::multifit::PosteriorOptimizerObjective)

// ignore std::vector methods that require a default ctor, as we don't have one for this class
%ignore std::vector<lsst::meas::multifit::PosteriorOptimizerIterationData>::vector(size_type);
%ignore std::vector<lsst::meas::multifit::PosteriorOptimizerIterationData>::resize(size_type);
%template(PosteriorOptimizerIterationDataVector) std::vector<lsst::meas::multifit::PosteriorOptimizerIterationData>;

%extend lsst::meas::multifit::PosteriorOptimizer {
    %pythoncode %{
        IterationData = PosteriorOptimizerIterationData
        Objective = PosteriorOptimizerObjective
        Control = PosteriorOptimizerControl

        def getConfig(self):
            config = self.ConfigClass()
            config.readControl(self.getControl())
            return config
    %}
}

%pythoncode %{
    import lsst.pex.config
    @lsst.pex.config.wrap(PosteriorOptimizerControl)
    class PosteriorOptimizerConfig(lsst.pex.config.Config):
        pass
    PosteriorOptimizer.ConfigClass = PosteriorOptimizerConfig
    PosteriorOptimizer.IterationDataVector = PosteriorOptimizerIterationDataVector
%}
