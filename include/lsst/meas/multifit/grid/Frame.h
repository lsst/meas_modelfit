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

#ifndef LSST_MEAS_MULTIFIT_GRID_Frame
#define LSST_MEAS_MULTIFIT_GRID_Frame

#include "lsst/meas/multifit/definition/Frame.h"
#include "lsst/meas/multifit/grid/Source.h"

namespace lsst { namespace meas { namespace multifit { namespace grid {

class Frame : public detail::FrameBase, private boost::noncopyable {
public:

    int const getPixelOffset() const { return _pixelOffset; }
    int const getPixelCount() const { return _pixelCount; }

    int const getFilterIndex() const { return _filterIndex; }
    int const getFrameIndex() const { return _frameIndex; }

    void applyWeights(ndarray::Array<double,2,1> const & matrix) const;

    void applyWeights(ndarray::Array<double,1,0> const & vector) const;

private:

    friend class Initializer;

    Frame(definition::Frame const & def, int pixelOffset, int filterIndex, int frameIndex);

    int _pixelOffset;
    int _pixelCount;

    int _filterIndex;
    int _frameIndex;
};

}}}} // namespace lsst::meas::multifit::grid

#endif // !LSST_MEAS_MULTIFIT_GRID_Frame
