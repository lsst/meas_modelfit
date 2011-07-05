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

//----------------------------- SharedElements -----------------------------------

%{
#include "lsst/meas/multifit/definition/SharedElement.h"
%}

SWIG_SHARED_PTR(grid_PositionElementPtr, lsst::meas::multifit::grid::SharedElement<lsst::meas::multifit::POSITION>);
SWIG_SHARED_PTR(grid_RadiusElementPtr, lsst::meas::multifit::grid::SharedElement<lsst::meas::multifit::RADIUS>);
SWIG_SHARED_PTR(grid_EllipticityElementPtr, lsst::meas::multifit::grid::SharedElement<lsst::meas::multifit::ELLIPTICITY>);

%define %DeclareGridSharedElement(TITLE, LOWER, UPPER, CONSTRAINT)
%template(grid_##TITLE##Element)
lsst::meas::multifit::grid::SharedElement<lsst::meas::multifit::UPPER>;
%PointerEQ(lsst::meas::multifit::grid::SharedElement<lsst::meas::multifit::UPPER>)
%AddStreamRepr(lsst::meas::multifit::grid::SharedElement<lsst::meas::multifit::UPPER>)
%extend lsst::meas::multifit::grid::SharedElement<lsst::meas::multifit::UPPER> {
    lsst::meas::multifit::TITLE const getValue() {
        return self->getValue();
    }
    bool isActive() {
        return self->isActive();
    }
}
%extend lsst::meas::multifit::grid::ObjectComponent {
    boost::shared_ptr< lsst::meas::multifit::grid::SharedElement< lsst::meas::multifit::UPPER > > get ## TITLE() {
        return self->get##TITLE();
    }
}
%enddef

%include "lsst/meas/multifit/grid/SharedElement.h"

//------------------------------------- ObjectComponent ---------------------------------------

%{
#include "lsst/meas/multifit/grid/ObjectComponent.h"
%}

%immutable lsst::meas::multifit::grid::ObjectComponent::sources;

SWIG_SHARED_PTR_DERIVED(grid_ObjectComponentPtr, lsst::meas::multifit::detail::ObjectComponentBase, lsst::meas::multifit::grid::ObjectComponent);

%rename(grid_ObjectComponent) lsst::meas::multifit::grid::ObjectComponent;
%immutable lsst::meas::multifit::grid::ObjectComponent::sources;

%include "lsst/meas/multifit/grid/ObjectComponent.h"

%DeclareGridSharedElement(Position, position, POSITION, detail::CircleConstraint);
%DeclareGridSharedElement(Radius, radius, RADIUS, detail::MinMaxConstraint);
%DeclareGridSharedElement(Ellipticity, ellipticity, ELLIPTICITY, detail::CircleConstraint);

%PointerEQ(lsst::meas::multifit::grid::ObjectComponent)
%AddStreamRepr(lsst::meas::multifit::grid::ObjectComponent)

%extend lsst::meas::multifit::grid::ObjectComponent {

    boost::shared_ptr< lsst::meas::multifit::grid::FluxGroup > getFluxGroup() {
        return self->getFluxGroup();
    }

}

//------------------------------------- Frame ---------------------------------------

%{
#include "lsst/meas/multifit/grid/Frame.h"
%}

SWIG_SHARED_PTR_DERIVED(grid_FramePtr, lsst::meas::multifit::detail::FrameBase, lsst::meas::multifit::grid::Frame);

%rename(grid_Frame) lsst::meas::multifit::grid::Frame;

%include "lsst/meas/multifit/grid/Frame.h"

%PointerEQ(lsst::meas::multifit::grid::Frame)
%AddStreamRepr(lsst::meas::multifit::grid::Frame)

%extend lsst::meas::multifit::grid::Frame {
    lsst::meas::multifit::FilterId getFilterId() { return self->getFilterId(); }
    lsst::meas::multifit::Wcs::Ptr getWcs() {
        return boost::const_pointer_cast<lsst::meas::multifit::Wcs>(self->getWcs());
    }
    lsst::meas::multifit::Psf::Ptr getPsf() {
        return boost::const_pointer_cast<lsst::meas::multifit::Psf>(self->getPsf());
    }
    lsst::meas::multifit::Footprint::Ptr getFootprint() {
        return boost::const_pointer_cast<lsst::meas::multifit::Footprint>(self->getFootprint());
    }
    lsst::ndarray::Array<lsst::meas::multifit::Pixel const,1,1> getData() { return self->getData(); }
    lsst::ndarray::Array<lsst::meas::multifit::Pixel const,1,1> getWeights() { return self->getWeights(); }
}

//------------------------------------- SourceComponent ---------------------------------------

%{
#include "lsst/meas/multifit/grid/SourceComponent.h"
%}

SWIG_SHARED_PTR(grid_SourceComponentPtr, lsst::meas::multifit::grid::SourceComponent);

%rename(grid_SourceComponent) lsst::meas::multifit::grid::SourceComponent;

%include "lsst/meas/multifit/grid/SourceComponent.h"

%PointerEQ(lsst::meas::multifit::grid::SourceComponent)

//------------------------------------- FluxGroup ---------------------------------------

%{
#include "lsst/meas/multifit/grid/FluxGroup.h"
%}

SWIG_SHARED_PTR_DERIVED(grid_FluxGroupPtr, lsst::meas::multifit::detail::FluxGroupBase, lsst::meas::multifit::grid::FluxGroup);
%rename(grid_FluxGroup) lsst::meas::multifit::grid::FluxGroup;
%immutable lsst::meas::multifit::grid::FluxGroup::components;

%include "lsst/meas/multifit/grid/FluxGroup.h"

//----------------------------- Grid -----------------------------------

%DeclareArray(lsst::meas::multifit::grid::SourceComponent, grid_SourceComponentArray, NO_INDEX)
%DeclareArray(lsst::meas::multifit::grid::ObjectComponent, grid_ObjectComponentArray, UNSORTED)
%DeclareArray(lsst::meas::multifit::grid::Frame, grid_FrameArray, SORTED)

%DeclareArray(lsst::meas::multifit::grid::PositionElement, grid_PositionArray, NO_INDEX)
%DeclareArray(lsst::meas::multifit::grid::RadiusElement, grid_RadiusArray, NO_INDEX)
%DeclareArray(lsst::meas::multifit::grid::EllipticityElement, grid_EllipticityArray, NO_INDEX)
%DeclareArray(lsst::meas::multifit::grid::FluxGroup, grid_FluxGroupArray, UNSORTED)

%DeclareArrayView(lsst::meas::multifit::grid::ObjectComponent, grid_FluxGroup_ComponentArray, SORTED)
%DeclareArrayView(lsst::meas::multifit::grid::SourceComponent, grid_ObjectComponent_SourceArray, NO_INDEX)

%{
#include "lsst/meas/multifit/grid/Grid.h"
%}

%immutable lsst::meas::multifit::grid::Grid::objects;
%immutable lsst::meas::multifit::grid::Grid::frames;
%immutable lsst::meas::multifit::grid::Grid::sources;
%immutable lsst::meas::multifit::grid::Grid::groups;
%immutable lsst::meas::multifit::grid::Grid::positions;
%immutable lsst::meas::multifit::grid::Grid::radii;
%immutable lsst::meas::multifit::grid::Grid::ellipticities;

SWIG_SHARED_PTR(GridPtr, lsst::meas::multifit::grid::Grid);

%include "lsst/meas/multifit/grid/Grid.h"
