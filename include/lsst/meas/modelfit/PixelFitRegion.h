// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2015-2016 LSST/AURA
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

#ifndef LSST_MEAS_MODELFIT_PixelFitRegion_h_INCLUDED
#define LSST_MEAS_MODELFIT_PixelFitRegion_h_INCLUDED

#include "lsst/pex/config.h"
#include "lsst/meas/modelfit/common.h"
#include "lsst/geom/Point.h"
#include "lsst/afw/image/Mask.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/geom/ellipses.h"

namespace lsst { namespace meas { namespace modelfit {


struct PixelFitRegionControl {

    PixelFitRegionControl() :
        nKronRadii(1.5),
        nPsfSigmaMin(4.0),
        nPsfSigmaGrow(2.0),
        nFitRadiiMin(1.0),
        nFitRadiiMax(3.0),
        maxArea(100000),
        maxBadPixelFraction(0.1)
    {
        badMaskPlanes.push_back("EDGE");
        badMaskPlanes.push_back("SAT");
        badMaskPlanes.push_back("BAD");
        badMaskPlanes.push_back("NO_DATA");
    }

    LSST_CONTROL_FIELD(
        nKronRadii, double,
        "Use this multiple of the Kron ellipse to set the fit region (for the final fit region, "
        "subject to the nFitRadiiMin and nFitRadiiMax constraints)."
    );

    LSST_CONTROL_FIELD(
        nPsfSigmaMin, double,
        "If the Kron radius is less than this multiple of the PSF width, ignore it and fall back "
        "to a PSF-oriented ellipse scaled to match the area of the footprint or this radius "
        "(whichever is larger)."
    );

    LSST_CONTROL_FIELD(
        nPsfSigmaGrow, double,
        "Grow the initial fit ellipses by this factor before comparing with the Kron/Footprint region"
    );

    LSST_CONTROL_FIELD(
        nFitRadiiMin, double,
        "Use this multiple of the initial fit ellipse then grow by the PSF width "
        "to determine the minimum final fit region size."
    );

    LSST_CONTROL_FIELD(
        nFitRadiiMax, double,
        "Use this multiple of the initial fit ellipse then grow by the PSF width "
        "to determine the maximum final fit region size."
    );

    LSST_CONTROL_FIELD(
        maxArea, int,
        "Abort if the fit region grows beyond this many pixels."
    );

    LSST_CONTROL_FIELD(
        badMaskPlanes, std::vector<std::string>,
        "Mask planes that indicate pixels that should be ignored in the fit."
    );

    LSST_CONTROL_FIELD(
        maxBadPixelFraction, double,
        "Maximum fraction of pixels that may be ignored due to masks; "
        "more than this and we don't even try."
    );

};


class PixelFitRegion {
public:

    PixelFitRegion(
        PixelFitRegionControl const & ctrl,
        afw::geom::ellipses::Quadrupole const & moments,
        afw::geom::ellipses::Quadrupole const & psfMoments,
        Scalar kronRadius,
        int footprintArea
    );

    PixelFitRegion(
        PixelFitRegionControl const & ctrl,
        afw::geom::ellipses::Quadrupole const & ellipse
    );

    bool applyEllipse(
        afw::geom::ellipses::Quadrupole const & deconvolved,
        afw::geom::ellipses::Quadrupole const & psfMoments
    );

    void applyMask(afw::image::Mask<> const & mask, geom::Point2D const & center);

    afw::geom::ellipses::Quadrupole ellipse;
    std::shared_ptr<afw::detection::Footprint> footprint;
    bool usedFootprintArea;
    bool usedPsfArea;
    bool maxArea;
    bool maxBadPixelFraction;
    bool usedMinEllipse;
    bool usedMaxEllipse;

private:
    PixelFitRegionControl _ctrl;
    afw::image::MaskPixel _badPixelMask;
};


}}} // lsst::meas::modelfit

#endif // !LSST_MEAS_MODELFIT_PixelFitRegion_h_INCLUDED
