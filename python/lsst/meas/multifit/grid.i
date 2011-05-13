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

//----------------------------- ParameterComponents -----------------------------------

%{
#include "lsst/meas/multifit/definition/parameters.h"
%}

SWIG_SHARED_PTR(grid_PositionComponentPtr, lsst::meas::multifit::grid::ParameterComponent<lsst::meas::multifit::POSITION>);
SWIG_SHARED_PTR(grid_RadiusComponentPtr, lsst::meas::multifit::grid::ParameterComponent<lsst::meas::multifit::RADIUS>);
SWIG_SHARED_PTR(grid_EllipticityComponentPtr, lsst::meas::multifit::grid::ParameterComponent<lsst::meas::multifit::ELLIPTICITY>);

%define %AddComponentAccessors(TITLE, LOWER, UPPER)
%enddef

%define %DeclareGridParameterComponent(TITLE, LOWER, UPPER, CONSTRAINT)
%template(grid_##TITLE##Component)
lsst::meas::multifit::grid::ParameterComponent<lsst::meas::multifit::UPPER>;
%PointerEQ(lsst::meas::multifit::grid::ParameterComponent<lsst::meas::multifit::UPPER>)
%AddStreamRepr(lsst::meas::multifit::grid::ParameterComponent<lsst::meas::multifit::UPPER>)
%extend lsst::meas::multifit::grid::ParameterComponent<lsst::meas::multifit::UPPER> {
    lsst::meas::multifit::TITLE const getValue() {
        return self->getValue();
    }
    bool isActive() {
        return self->isActive();
    }
}
%extend lsst::meas::multifit::grid::Object {
    boost::shared_ptr< lsst::meas::multifit::grid::ParameterComponent< lsst::meas::multifit::UPPER > > get ## TITLE() {
        return self->get##TITLE();
    }
}
%enddef

%include "lsst/meas/multifit/grid/parameters.h"

//-------------------------------------- Array -----------------------------------------

%{
#include "lsst/meas/multifit/grid/Array.h"
%}

namespace lsst { namespace meas { namespace multifit { namespace grid {

template <typename T> class Array {};

}}}} // namespace lsst::meas::multifit::grid

%define %DeclareArray(NAME)
%template(grid_##NAME##Array)
lsst::meas::multifit::grid::Array<lsst::meas::multifit::grid::NAME>;
%extend lsst::meas::multifit::grid::Array<lsst::meas::multifit::grid::NAME> {
    NAME const * __getitem__(int n) {
        if (n < 0 || n >= self->size()) {
            PyErr_SetString(PyExc_IndexError, "Index out of range.");
            return 0;
        }
        return &self->operator[](n);
    }
    int __len__() {
        return self->size();
    }
    %pythoncode %{
        def __iter__(self):
            for id in range(len(self)):
                yield self[id]
        def __repr__(self):
            return "NAME""Array([\n"  + ",".join(repr(v) for v in self) + "\n])"
    %}
}
%enddef

//------------------------------------- Object ---------------------------------------

%{
#include "lsst/meas/multifit/grid/Object.h"
%}

%immutable lsst::meas::multifit::grid::Object::sources;

SWIG_SHARED_PTR_DERIVED(grid_ObjectPtr, lsst::meas::multifit::detail::ObjectBase, lsst::meas::multifit::grid::Object);

%rename(grid_Object) lsst::meas::multifit::grid::Object;

%include "lsst/meas/multifit/grid/Object.h"

%DeclareGridParameterComponent(Position, position, POSITION, detail::CircleConstraint);
%DeclareGridParameterComponent(Radius, radius, RADIUS, detail::MinMaxConstraint);
%DeclareGridParameterComponent(Ellipticity, ellipticity, ELLIPTICITY, detail::CircleConstraint);

%PointerEQ(lsst::meas::multifit::grid::Object)
%AddStreamRepr(lsst::meas::multifit::grid::Object)

%extend lsst::meas::multifit::grid::Object {

    lsst::afw::geom::Point2D makePoint(lsst::ndarray::Array<double const,1,1> const & params) const {
        return self->makePoint(params.getData());
    }

    lsst::afw::geom::ellipses::Ellipse makeEllipse(lsst::ndarray::Array<double const,1,1> const & params) const {
        return self->makeEllipse(params.getData());
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

//------------------------------------- Source ---------------------------------------

%{
#include "lsst/meas/multifit/grid/Source.h"
%}

SWIG_SHARED_PTR(grid_SourcePtr, lsst::meas::multifit::grid::Source);

%rename(grid_Source) lsst::meas::multifit::grid::Source;

%include "lsst/meas/multifit/grid/Source.h"

%PointerEQ(lsst::meas::multifit::grid::Source)

//----------------------------- Grid -----------------------------------

%DeclareArray(Source)
%DeclareArray(Object)
%DeclareArray(Frame)

%{
#include "lsst/meas/multifit/grid/Grid.h"
%}

%immutable lsst::meas::multifit::grid::Grid::objects;
%immutable lsst::meas::multifit::grid::Grid::frames;
%immutable lsst::meas::multifit::grid::Grid::sources;

SWIG_SHARED_PTR(GridPtr, lsst::meas::multifit::grid::Grid);

%include "lsst/meas/multifit/grid/Grid.h"
