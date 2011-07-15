// -*- LSST-C++ -*-
/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010, 2011 LSST Corporation.
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

#ifndef LSST_MEAS_MULTIFIT_GRID_FluxGroup
#define LSST_MEAS_MULTIFIT_GRID_FluxGroup

#include "lsst/meas/multifit/definition/FluxGroup.h"
#include "lsst/meas/multifit/grid/SourceComponent.h"
#include "lsst/meas/multifit/grid/ObjectComponent.h"
#include "lsst/meas/multifit/grid/Frame.h"
#include "lsst/meas/multifit/containers/Array.h"
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

namespace lsst { namespace meas { namespace multifit { namespace grid {

class FluxGroup : public detail::FluxGroupBase {
public:
    
    typedef boost::shared_ptr< FluxGroup > Ptr;
    typedef boost::shared_ptr< FluxGroup const > ConstPtr;

    typedef containers::ArrayView<ObjectComponent,containers::SORTED> ComponentArray;

    /// @brief Array of ObjectComponents in this group.
    ComponentArray components;

    /// @brief Number of coefficients in the flux group per instance (frame/filter).
    int const getSourceCoefficientCount() const { return _sourceCoefficientCount; }

    /// @brief Return the coefficient offset for the frame/filter with the given index.
    int const getCoefficientOffset(int frameOrFilterIndex) const {
        return _coefficientOffset + _sourceCoefficientCount * frameOrFilterIndex;
    }

    /// @brief Return the coefficient offset for the given frame.
    int const getCoefficientOffset(Frame const & frame) const {
        return getCoefficientOffset(isVariable() ? frame.getFrameIndex() : frame.getFilterIndex());
    }

    /// @brief Return the constraint offset for the frame/filter with the given index.
    int const getConstraintOffset(int frameOrFilterIndex) const {
        return _constraintOffset + getConstraintCount() * frameOrFilterIndex;
    }

    /// @brief Return the coefficient offset for the given frame.
    int const getConstraintOffset(Frame const & frame) const {
        return getConstraintOffset(isVariable() ? frame.getFrameIndex() : frame.getFilterIndex());
    }

    //@{
    /// @brief Return linear inequality constraints for one instance (frame/filter) of this group.
    lsst::ndarray::Array<Pixel const,2,2> getConstraintMatrix() const { return _constraintMatrix; }
    lsst::ndarray::Array<Pixel const,1,1> getConstraintVector() const { return _constraintVector; }
    int const getConstraintCount() const { return _constraintVector.getSize<0>(); }
    //@}

    /// @brief Return an array that computes the flux for one instance (frame/filter) of this group.
    lsst::ndarray::Array<Pixel const,1,1> getIntegration() const { return _integration; }

#ifndef SWIG // these are wrapped explicitly; SWIG is confused by the typedefs and "bool &"
    // Use const accessors from base class.
    using detail::FluxGroupBase::isVariable;
    using detail::FluxGroupBase::getMaxMorphologyRatio;
#endif 

private:
    
    friend class Initializer;

    FluxGroup(definition::FluxGroup const & definition, int coefficientOffset) :
        detail::FluxGroupBase(definition),
        _coefficientOffset(coefficientOffset),
        _sourceCoefficientCount(0), _constraintOffset(0)
    {}

    // Sets up constraints; must be called after components array view is initialized.
    void initialize();

    int _coefficientOffset;
    int _sourceCoefficientCount;
    int _constraintOffset;
    ndarray::Array<Pixel,2,2> _constraintMatrix;
    ndarray::Array<Pixel,1,1> _constraintVector;
    ndarray::Array<Pixel,1,1> _integration;
};

#ifndef SWIG
inline afw::geom::Point2D const SourceComponent::getReferencePoint() const {
    return _transform(object.getPosition()->getValue());
}

inline int const SourceComponent::getCoefficientOffset() const {
    return object.getFluxGroup()->getCoefficientOffset(frame) + object.getGroupCoefficientOffset();
}

inline int const SourceComponent::getCoefficientCount() const {
    return object.getSourceCoefficientCount();
}
#endif

}}}} // namespace lsst::meas::multifit::grid

#endif // !LSST_MEAS_MULTIFIT_GRID_FluxGroup
