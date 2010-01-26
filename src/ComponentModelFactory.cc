// -*- lsst-c++ -*-
#include "lsst/meas/multifit/ComponentModel.h"
#include "lsst/meas/multifit/ComponentModelFactory.h"

namespace multifit = lsst::meas::multifit;

// ComponentModelFactory -------------------------------------------------------
multifit::Model::Ptr multifit::ComponentModelFactory::makeModel(
    ParameterVector const & linearParameters,
    ParameterVector const & nonlinearParameters
) const {
    Model::Ptr model(
        new ComponentModel(
            linearParameters.size(), _astrometryTemplate, _morphologyTemplate
        )
    );
    model->setLinearParameters(linearParameters);
    model->setNonlinearParameters(nonlinearParameters);
    return model;
}
