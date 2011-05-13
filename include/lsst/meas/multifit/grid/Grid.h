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
#include "lsst/meas/multifit/definition/Definition.h"

#include <boost/scoped_array.hpp>
#include <map>

namespace lsst { namespace meas { namespace multifit { namespace grid {

class Grid {
public:
    typedef boost::shared_ptr<Grid> Ptr;

    typedef grid::Object Object;
    typedef grid::Frame Frame;
    typedef grid::Source Source;

    typedef grid::Array<Object> ObjectArray;
    typedef grid::Array<Frame> FrameArray;
    typedef grid::Array<Source> SourceArray;

    typedef grid::ComponentArray<POSITION> PositionArray;
    typedef grid::ComponentArray<RADIUS> RadiusArray;
    typedef grid::ComponentArray<ELLIPTICITY> EllipticityArray;

    Definition makeDefinition() const;

    Definition makeDefinition(double const * paramIter) const;

#ifndef SWIG
    void writeParameters(double * paramIter) const;
#endif
    void writeParameters(lsst::ndarray::Array<double, 1, 1> const & params) const {
        writeParameters(params.getData());
    }

    double sumLogWeights() const;

    int const getFilterIndex(FilterId filterId) const;

    /// @brief Return true if all parameters are in-bounds.
    bool checkBounds(double const * paramIter) const;

    /**
     *  @brief Clip any out-of-bounds parameters to the bounds and return a positive number
     *         indicating how much clipping was necessary.
     *
     *  The returned value has no well-defined units and may penalize some parameter types
     *  more than others.  The return value will be zero when no clipping is necessary.
     */
    double clipToBounds(double * paramIter) const;

    ~Grid();

    ObjectArray objects;
    FrameArray frames;
    SourceArray sources;

#ifndef SWIG
    //@{
    /// Arrays of all active parameter components (inactive ones are still held by the Objects). 
    PositionArray positions;
    RadiusArray radii;
    EllipticityArray ellipticities;
    //@}
#endif

    int const getFilterCount() const { return _filterCount; }
    int const getCoefficientCount() const { return _coefficientCount; }
    int const getPixelCount() const { return _pixelCount; }
    int const getParameterCount() const { return _parameterCount; }

    CONST_PTR(Wcs) const getWcs() const { return _wcs; }

    static Ptr make(Definition const & definition) {
        return Ptr(new Grid(definition));
    }

private:

    friend class Initializer;

    typedef std::map<FilterId, int > FilterMap;

    explicit Grid(Definition const & definition);

    void _destroy();

    template <typename ObjectIterator, typename FrameIterator>
    void _initialize(
        ObjectIterator const & objectBegin, ObjectIterator const & objectEnd,
        FrameIterator const & frameBegin, FrameIterator const & frameEnd
    );

    int _filterCount;
    int _coefficientCount;
    int _pixelCount;
    int _parameterCount;

    boost::scoped_array<char> _objectData;
    boost::scoped_array<char> _frameData;
    boost::scoped_array<char> _sourceData;
    CONST_PTR(Wcs) _wcs;
    FilterMap _filters;
};

} // namespace grid

typedef grid::Grid Grid;

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_Grid
