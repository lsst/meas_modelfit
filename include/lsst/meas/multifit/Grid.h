// -*- LSST-C++ -*-
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

#ifndef LSST_MEAS_MULTIFIT_Grid
#define LSST_MEAS_MULTIFIT_Grid

#include "lsst/meas/multifit/grid/Frame.h"
#include "lsst/meas/multifit/grid/Object.h"
#include "lsst/meas/multifit/grid/Array.h"
#include "lsst/meas/multifit/Definition.h"

#include <boost/scoped_array.hpp>
#include <map>

namespace lsst { namespace meas { namespace multifit {

class Grid {
public:

    typedef grid::Array<grid::Object> ObjectArray;
    typedef grid::Array<grid::Frame> FrameArray;
    typedef grid::Array<grid::Source> SourceArray;

    typedef grid::ComponentArray<grid::PositionComponent> PositionArray;
    typedef grid::ComponentArray<grid::RadiusComponent> RadiusArray;
    typedef grid::ComponentArray<grid::EllipticityComponent> EllipticityArray;

    typedef std::map< definition::Filter const *, int > FilterMap;

    explicit Grid(Definition const & definition);

    Grid(Grid const & other);

    Definition makeDefinition() const;
    Definition makeDefinition(double const * param_iter) const;

    void writeParameters(double * param_iter) const;

    ~Grid() { _destroy(); }

    ObjectArray objects;
    FrameArray frames;
    SourceArray sources;

    PositionArray positions;
    RadiusArray radii;
    EllipticityArray ellipticities;

    FilterMap filters;

    int filter_count;
    int coefficient_count;
    int pixel_count;
    int parameter_count;

    Wcs::Ptr wcs;

private:

    void _destroy();

    template <typename ObjectIterator, typename FrameIterator>
    void _initialize(
        ObjectIterator const & object_begin, ObjectIterator const & object_end,
        FrameIterator const & frame_begin, FrameIterator const & frame_end
    );

    boost::scoped_array<char> _object_data;
    boost::scoped_array<char> _frame_data;
    boost::scoped_array<char> _source_data;

};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_Grid
