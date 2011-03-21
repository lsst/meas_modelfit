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

#ifndef LSST_MEAS_MULTIFIT_constants
#define LSST_MEAS_MULTIFIT_constants

#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/image/Wcs.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/image/Calib.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/LocalPsf.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/detection/Source.h"
#include "lsst/pex/exceptions.h"

#ifndef SWIG

#define FRIEND_MAKE_SHARED_1(T, A1)                                       \
    friend boost::shared_ptr<T> boost::make_shared<T,A1>(A1 const &)

#define FRIEND_MAKE_SHARED_2(T, A1, A2)                                   \
    friend boost::shared_ptr<T> boost::make_shared<T,A1,A2>(A1 const &, A2 const &)

#define FRIEND_MAKE_SHARED_3(T, A1, A2, A3)                               \
    friend boost::shared_ptr<T> boost::make_shared<T,A1,A2,A3>(A1 const &, A2 const &, A3 const &)

#define FRIEND_MAKE_SHARED_4(T, A1, A2, A3, A4)                         \
    friend boost::shared_ptr<T> boost::make_shared<T,A1,A2,A3,A4>(      \
        A1 const &, A2 const &, A3 const &, A4 const &                  \
    )
#endif

namespace lsst { namespace meas { namespace multifit {

namespace detail {

inline void checkSize(int actualSize, int expectedSize, char const * message) {
    if (actualSize != expectedSize) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LengthErrorException, 
            (boost::format(message) % actualSize % expectedSize).str()
        );
    }
}

} // namespace detail



typedef boost::int64_t ID;
typedef int FilterId;

enum ParameterType { POSITION, RADIUS, ELLIPTICITY };

typedef lsst::afw::geom::ellipses::TraceRadius Radius;
typedef lsst::afw::geom::ellipses::LogShear Ellipticity;
typedef lsst::afw::geom::ellipses::Separable<Ellipticity, Radius> EllipseCore;
typedef lsst::afw::geom::ellipses::Ellipse Ellipse;

typedef lsst::afw::image::Filter Filter;
typedef lsst::afw::image::Calib Calib;
typedef lsst::afw::image::Wcs Wcs;
typedef lsst::afw::detection::Psf Psf;
typedef lsst::afw::detection::LocalPsf LocalPsf;
typedef lsst::afw::detection::Footprint Footprint;

typedef LocalPsf::Pixel Pixel;
typedef lsst::afw::image::Exposure<lsst::meas::multifit::Pixel> Exposure;

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_constants
