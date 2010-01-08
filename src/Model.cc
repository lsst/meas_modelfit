#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/ModelProjection.h"

namespace multifit = lsst::meas::multifit;

void multifit::Model::setLinearParameters(ParameterConstIterator const parameters) {
    std::copy(parameters, parameters + getLinearParameterSize(), _linearParameters->data());
    _handleLinearParameterChange();
    _broadcastLinearParameterChange();
}

void multifit::Model::setNonlinearParameters(ParameterConstIterator const parameters) {
    std::copy(parameters, parameters + getNonlinearParameterSize(), _nonlinearParameters->data());
    _handleNonlinearParameterChange();
    _broadcastNonlinearParameterChange();
}

void multifit::Model::_broadcastLinearParameterChange() const {
    ProjectionList::iterator const & end(_projectionList.end());
    ProjectionList::iterator i(_projectionList.begin());
    while (i != end) {
        ModelProjection::Ptr projection(i->lock());
        if (!projection) {
            i = _projectionList.erase(i);
        } else {
            projection->_handleLinearParameterChange();
            ++i;
        }
    }
}

void multifit::Model::_broadcastNonlinearParameterChange() const {
    ProjectionList::iterator const & end(_projectionList.end());
    ProjectionList::iterator i = _projectionList.begin();
    while (i != end) {
        ModelProjection::Ptr projection(i->lock());
        if (!projection) {
            i = _projectionList.erase(i);
        } else {
            projection->_handleNonlinearParameterChange();
            ++i;
        }
    }
}

void multifit::Model::_registerProjection(
    boost::shared_ptr<ModelProjection> const & projection
) const {
    ProjectionWeakPtr weak(projection);
    _projectionList.push_back(weak);
}
