#include "lsst/meas/multifit/components/MorphologyProjection.h"
#include "lsst/meas/multifit/components/Morphology.h"

namespace components = lsst::meas::multifit::components;

int const components::MorphologyProjection::getLinearParameterSize() const {
    return _morphology->getLinearParameterSize();
}

int const components::MorphologyProjection::getNonlinearParameterSize() const {
    return _morphology->getNonlinearParameterSize();
}
