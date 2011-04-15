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
    typedef boost::shared_ptr<Grid> Ptr;
    typedef boost::shared_ptr<const Grid> ConstPtr;

    typedef grid::Array<grid::Object> ObjectArray;
    typedef grid::Array<grid::Frame> FrameArray;
    typedef grid::Array<grid::Source> SourceArray;

    typedef grid::ComponentArray<grid::PositionComponent> PositionArray;
    typedef grid::ComponentArray<grid::RadiusComponent> RadiusArray;
    typedef grid::ComponentArray<grid::EllipticityComponent> EllipticityArray;

    typedef std::map<FilterId, int > FilterMap;

    explicit Grid(Definition const & definition);

    Grid(Grid const & other);

    Definition makeDefinition() const;
    Definition makeDefinition(double const * paramIter) const;

    void writeParameters(double * paramIter) const;

    double sumLogWeights() const;

    ~Grid() { _destroy(); }

    ObjectArray objects;
    FrameArray frames;
    SourceArray sources;

    PositionArray positions;
    RadiusArray radii;
    EllipticityArray ellipticities;

    FilterMap filters;
    
    int filterCount;
    int coefficientCount;
    int pixelCount;
    int parameterCount;

    Wcs::Ptr wcs;

private:

    void _destroy();

    template <typename ObjectIterator, typename FrameIterator>
    void _initialize(
        ObjectIterator const & objectBegin, ObjectIterator const & objectEnd,
        FrameIterator const & frameBegin, FrameIterator const & frameEnd
    );

    boost::scoped_array<char> _objectData;
    boost::scoped_array<char> _frameData;
    boost::scoped_array<char> _sourceData;

};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_Grid
