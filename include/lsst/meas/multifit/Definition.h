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

#ifndef LSST_MEAS_MULTIFIT_Definition
#define LSST_MEAS_MULTIFIT_Definition

#include "lsst/meas/multifit/definition/Frame.h"
#include "lsst/meas/multifit/definition/Object.h"
#include "lsst/meas/multifit/definition/Set.h"

namespace lsst { namespace meas { namespace multifit {

class Definition {
public:
    typedef definition::Frame Frame;
    typedef definition::Object Object;

    typedef definition::Set<Frame,ID,&Frame::id> FrameSet;
    typedef definition::Set<Object,ID,&Object::id> ObjectSet; 

    Definition() {}

    explicit Definition(Wcs::Ptr const & wcs_) : wcs(wcs_) {}

    Definition(Definition const & other);

    FrameSet frames;
    ObjectSet objects;

    Wcs::Ptr wcs;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_Definition
