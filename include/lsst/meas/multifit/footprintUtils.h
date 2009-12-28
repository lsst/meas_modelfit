#ifndef LSST_MEAS_MULTIFIT_FOOTPRINT_UTILS_H
#define LSST_MEAS_MULTIFIT_FOOTPRINT_UTILS_H

#include "boost/shared_ptr.hpp"

#include "lsst/afw/geom/Box.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/geom/ellipses.h"

lsst::afw::detection::Footprint::Ptr makeFootprint(
    lsst::afw::geom::ellipses::Ellipse const & ellipse
) {
    lsst::afw::detection::Footprint::Ptr fp(new lsst::afw::detection::Footprint());
    lsst::afw::geom::ellipses::Core::RadialFraction rf(ellipse.getCore());

    lsst::afw::geom::BoxD envelope = ellipse.computeEnvelope();
    int const yEnd = envelope.getMaxY() + 1;
    int const xEnd = envelope.getMaxX() + 1;
    lsst::afw::geom::ExtentD dp(lsst::afw::geom::PointD(0) -ellipse.getCenter());
    for (int y = envelope.getMinY(); y<yEnd; ++y) {
        int x = envelope.getMinX();
        while (rf(lsst::afw::geom::PointD::makeXY(x,y) + dp) > 1.0) {
            if (x >= xEnd) {
                if (++y >= yEnd) 
                    return fp;
                x = envelope.getMinX();
            } else {
                ++x;
            }
        }
        int start = x;
        while (rf(lsst::afw::geom::PointD::makeXY(x,y) + dp) <= 1.0 && x < xEnd) 
            ++x;

        fp->addSpan(y, start, x-1);
    }
    return fp;
}

#endif
