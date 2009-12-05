#include <ndarray.hpp>
#include <lsst/meas/multifit/Model.h>
#include <lsst/meas/multifit/ModelFactory.h>
#include <lsst/meas/multifit/projections/ModelProjection.h>

namespace multifit = lsst::meas::multifit;

multifit::Model(int linearParameterSize, int nonlinearParameterSize) :
    _linearParameterVector(new ParameterVector(linearParameterSize)),
    _nonlinearParameterVector(new ParameterVector(nonlinearParameterSize)),
    _projectionList()
{}

multifit::Model(Model const & model) :
    _linearParameterVector(new ParameterVector(model.getLinearParameters())),
    _nonlinearParameterVector(
        new ParameterVector(model.getNonlinearParameters())
    ),
    _projectionList()
{}

multifit::projections::ModelProjection::Ptr multifit::Model::makeProjection(
    boost::shared_ptr<const Kernel> const & kernel,
    boost::shared_ptr<const Wcs> const & wcs,
    boost::shared_ptr<const Footprint> const & footprint,
    double photFactor,
    int activeProducts 
) const {
    projections::ModelProjection::Ptr projection(createProjection());
    projection->_model = shared_from_this();
    projection->_wcs = wcs;
    projection->_footprint = footprint;
    projection->_photFactor = photFactor;
    projection->setKernel(kernel);
    projection->enableProducts(activeProducts);
    projection->_handleLinearParameterChange();
    projection->_handleNonlinearParameterChange();
    _projectionList.push_back(projection);
    return projection;
}

void multifit::Model::setLinearParameters(
    ParameterConstIterator const parameters
) {
    std::copy(
        parameters,
        parameters + getLinearParameterSize(),
        _linearParameterVector->data()
    );
    _handleLinearParameterChange();
    _broadcastLinearParameterChange();
}


void multifit::Model::setNonlinearParameters(
    ParameterConstIterator const parameters
) {
    std::copy(
        parameters,
        parameters + getNonlinearParameterSize(),
        _nonlinearParameterVector->data()
    );
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
//---------------------------------------------------------------------------//
