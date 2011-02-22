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

#ifndef LSST_MEAS_MULTIFIT_GRID_Source
#define LSST_MEAS_MULTIFIT_GRID_Source

#include "lsst/afw/geom/AffineTransform.h"
#include "lsst/meas/multifit/ModelBasis.h"

namespace lsst { namespace meas { namespace multifit { namespace grid {

class Frame;
class Object;

class Source {
public:

    Source(
        Frame const & frame, Object const & object, 
        CONST_PTR(afw::image::Wcs) const & wcs
    );

    Frame const & frame;
    Object const & object;

    afw::geom::AffineTransform transform;
    ModelBasis::Ptr basis;
    afw::detection::LocalPsf::Ptr localPsf;
    mutable void * extra;

};

}}}} // namespace lsst::meas::multifit::grid

#endif // !LSST_MEAS_MULTIFIT_GRID_sources
