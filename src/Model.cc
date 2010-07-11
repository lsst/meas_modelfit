// -*- lsst-c++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/ModelProjection.h"

namespace multifit = lsst::meas::multifit;

/**
 * @name Model Parameter Setters
 *
 * Set the parameters of this model. 
 * All registered projections will be notified of change
 */
//@{
void multifit::Model::setLinearParameters(
    ParameterConstIterator const parameters
) {
    std::copy(
        parameters, 
        parameters + getLinearParameterSize(), 
        _linearParameters->data()
    );
    _handleLinearParameterChange();
    _broadcastLinearParameterChange();
}
void multifit::Model::setLinearParameters(ParameterVector const & parameters) {
    setLinearParameters(parameters.data());
}

void multifit::Model::setNonlinearParameters(
    ParameterConstIterator const parameters
) {
    std::copy(
        parameters, 
        parameters + getNonlinearParameterSize(), 
        _nonlinearParameters->data()
    );
    _handleNonlinearParameterChange();
    _broadcastNonlinearParameterChange();
}
void multifit::Model::setNonlinearParameters(ParameterVector const & parameters) {
    setNonlinearParameters(parameters.data());
}
//@}

/**
 * Notify all associated ModelProjections that the linear parameters have changed.
 */
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

/**
 * Notify all associated ModelProjections that the nonlinear parameters have changed.
 */
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

/**
 * Add a projection to the list of registered listeners.
 *
 * All registered Model Projections will be notified when the linear or nonlinear
 * parameters of this Model are modified
 */
void multifit::Model::_registerProjection(
    boost::shared_ptr<ModelProjection> const & projection
) const {
    ProjectionWeakPtr weak(projection);
    _projectionList.push_back(weak);
}
