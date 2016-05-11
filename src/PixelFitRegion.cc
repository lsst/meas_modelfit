// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2015 LSST/AURA.
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
#include <cmath>

#include "lsst/meas/modelfit/PixelFitRegion.h"

namespace lsst { namespace meas { namespace modelfit {

namespace {

afw::image::MaskPixel initBadPixelMask(std::vector<std::string> const & planes) {
    afw::image::MaskPixel badPixelMask = 0x0;
    for (std::vector<std::string>::const_iterator iter = planes.begin(); iter != planes.end(); ++iter) {
        badPixelMask |= afw::image::Mask<>::getPlaneBitMask(*iter);
    }
    return badPixelMask;
}

} // anonymous


PixelFitRegion::PixelFitRegion(
    PixelFitRegionControl const & ctrl,
    afw::geom::ellipses::Quadrupole const & moments,
    afw::geom::ellipses::Quadrupole const & psfMoments,
    Scalar kronRadius,
    int footprintArea
) : ellipse(moments),
    usedFootprintArea(false),
    usedPsfArea(false),
    maxArea(false),
    maxBadPixelFraction(false),
    usedMinEllipse(false),
    usedMaxEllipse(false),
    _ctrl(ctrl),
    _badPixelMask(initBadPixelMask(ctrl.badMaskPlanes))
{
    // Try setting ellipse to a multiple of the Kron ellipse, fall back
    // to Footprint area circle and then scaled PSF moments if necessary.
    Scalar psfArea = _ctrl.nPsfSigmaMin*_ctrl.nPsfSigmaMin*psfMoments.getArea();
    try {
        ellipse.normalize();
        if (kronRadius > 0.0) {
            ellipse.scale(_ctrl.nKronRadii * kronRadius / ellipse.getDeterminantRadius());
            if (ellipse.getArea() < psfArea) {
                usedFootprintArea = true;
            }
        } else {
            usedFootprintArea = true;
        }
    } catch (pex::exceptions::InvalidParameterError &) {
        usedFootprintArea = true;
    }
    if (usedFootprintArea) {
        ellipse = psfMoments;
        if (psfArea > footprintArea) {
            // can't actually use footprint area after all, it's too small
            usedFootprintArea = false;
            usedPsfArea = true;
            ellipse.scale(_ctrl.nPsfSigmaMin);
        } else {
            ellipse.scale(std::sqrt(footprintArea/psfMoments.getArea()));
        }
    }
    if (ellipse.getArea() > _ctrl.maxArea) {
        maxArea = true;
    }
}


PixelFitRegion::PixelFitRegion(
    PixelFitRegionControl const & ctrl,
    afw::geom::ellipses::Quadrupole const & ellipse
) : ellipse(ellipse),
    usedFootprintArea(false),
    usedPsfArea(false),
    maxArea(false),
    maxBadPixelFraction(false),
    usedMinEllipse(false),
    usedMaxEllipse(false),
    _ctrl(ctrl),
    _badPixelMask(initBadPixelMask(ctrl.badMaskPlanes))
{
    if (ellipse.getArea() > _ctrl.maxArea) {
        maxArea = true;
    }
}


bool PixelFitRegion::applyEllipse(
    afw::geom::ellipses::Quadrupole const & deconvolved,
    afw::geom::ellipses::Quadrupole const & psfMoments
) {
    bool constrained = false;

    // Create a scaled PSF ellipse we'll use to convolve (i.e. add in quadrature).
    afw::geom::ellipses::Quadrupole ePsfGrow(psfMoments);
    ePsfGrow.scale(_ctrl.nPsfSigmaGrow);

    // alpha: factor to scale the deconvolved ellipse before convolving it with the PSF.
    Scalar alpha = 0.0;

    if (_ctrl.nFitRadiiMin == _ctrl.nFitRadiiMax) {
        // alpha is completely constrained by the fit radius
        alpha = _ctrl.nFitRadiiMin;
        usedMinEllipse = true;
        usedMaxEllipse = true;
        constrained = true;
    } else {
        // We solve an quadratic equation to find the alpha that gives an ellipse with the same
        // area as the current one (but with ellipticity determined by the given
        // deconvolved ellipse).
        Scalar a = deconvolved.getDeterminant();
        Scalar b = deconvolved.getIxx()*ePsfGrow.getIyy() + deconvolved.getIyy()*ePsfGrow.getIxx()
            - 2*deconvolved.getIxy()*ePsfGrow.getIxy();
        Scalar c = ePsfGrow.getDeterminant() - ellipse.getDeterminant();
        alpha = (std::sqrt(b*b - 4.0*a*c) - b) / (2.0*a);
        if (alpha < _ctrl.nFitRadiiMin) {
            alpha = _ctrl.nFitRadiiMin;
            usedMinEllipse = true;
            constrained = true;
        } else if (alpha > _ctrl.nFitRadiiMax) {
            alpha = _ctrl.nFitRadiiMax;
            usedMaxEllipse = true;
            constrained = true;
        }
    }

    // Regardless of how we set alpha, we always use the deconvolved ellipse to set the ellipticity,
    // then grow by the PSF size.
    ellipse = deconvolved;
    ellipse.scale(alpha);
    ellipse.convolve(ePsfGrow).inPlace();

    return constrained;
}

void PixelFitRegion::applyMask(afw::image::Mask<> const & mask, afw::geom::Point2D const & center) {
    Scalar originalArea = ellipse.getArea();
    footprint = std::make_shared<afw::detection::Footprint>(
        afw::geom::ellipses::Ellipse(ellipse, center),
        mask.getBBox(lsst::afw::image::PARENT)
    );
    footprint->clipTo(mask.getBBox(afw::image::PARENT));
    if (footprint->getArea() == 0) {
        maxBadPixelFraction = true;
        footprint.reset();
        return;
    }
    footprint->intersectMask(mask, _badPixelMask);
    if (originalArea - footprint->getArea() > originalArea*_ctrl.maxBadPixelFraction) {
        maxBadPixelFraction = true;
        footprint.reset();
        return;
    }
}




}}} // lsst::meas::modelfit
