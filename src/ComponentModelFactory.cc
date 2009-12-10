#include "lsst/meas/multifit/ComponentModel.h"
#include "lsst/meas/multifit/ComponentModelFactory.h"

namespace multifit = lsst::meas::multifit;

// ComponentModelFactory -------------------------------------------------------
multifit::Model::Ptr multifit::ComponentModelFactory::makeModel(
    int linearParameterSize,
    ParameterConstIterator linearParameterIter,
    ParameterConstIterator nonlinearParameterIter
) const {
    Model::Ptr model(
        new ComponentModel(
            linearParameterSize,
            _astrometryTemplate,
            _morphologyTemplate)
    );
    model->setLinearParameters(linearParameterIter);
    model->setNonlinearParameters(nonlinearParameterIter);
    return model;
}
