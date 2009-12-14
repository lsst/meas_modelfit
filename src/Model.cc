#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/ModelProjection.h"

namespace multifit = lsst::meas::multifit;

void multifit::Model::setLinearParameters(ParameterConstIterator const parameters) {
    std::copy(parameters, parameters + getLinearParameterSize(), _linearParameterVector->data());
    _handleLinearParameterChange();
    _broadcastLinearParameterChange();
}

void multifit::Model::setNonlinearParameters(ParameterConstIterator const parameters) {
    std::copy(parameters, parameters + getNonlinearParameterSize(), _nonlinearParameterVector->data());
    _handleNonlinearParameterChange();
    _broadcastNonlinearParameterChange();
}

void multifit::Model::broadcastLinearParameterChange() {
    ProjectionList::iterator const & end(_projectionList.end());
    ProjectionList::iterator i(_projectionList.begin());
    while (i != end) {
        projections::ModelProjection::Ptr projection(i->lock());
        if (!projection) {
            i = _projectionList.erase(i);
        } else {
            projection->_handleLinearParameterChange();
            ++i;
        }
    }
}

void multifit::Model::broadcastNonlinearParameterChange() {
    ProjectionList::iterator const & end(_projectionList.end());
    ProjectionList::iterator i = _projectionList.begin();
    while (i != end) {
        projections::ModelProjection::Ptr projection(i->lock());
        if (!projection) {
            i = _projectionList.erase(i);
        } else {
            projection->_handleNonlinearParameterChange();
            ++i;
        }
    }
}

void multifit::Model::_registerProjection(ModelProjection::Ptr const & projection) const {
    ProjectionWeakPtr weak(projection);
    _projectionList.push_back(weak);
}
