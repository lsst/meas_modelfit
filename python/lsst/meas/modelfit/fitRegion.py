#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2013 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import lsst.pex.config
import lsst.pipe.base
import lsst.afw.image
import lsst.afw.detection
import lsst.afw.geom.ellipses

__all__ = ("FitRegionConfig", "setupFitRegion")

class FitRegionConfig(lsst.pex.config.Config):
    """Configuration for which pixels to include in a fit.

    See setupFitRegion for more information.
    """
    nGrowFootprint = lsst.pex.config.Field(
        dtype=int,
        default=0,
        optional=True,
        doc="How many pixels to grow detection footprints by (may be negative to shrink)."
    )
    nInitialRadii = lsst.pex.config.Field(
        dtype=float,
        default=3.0,
        optional=True,
        doc="Scale the initial model ellipse by this much."
    )
    badMaskPlanes = lsst.pex.config.ListField(
        dtype=str,
        default=["SAT", "NO_DATA"],
        doc="Mask planes that indicate pixels that should not be included in the fit."
    )

    def validate(self):
        lsst.pex.config.Config.validate(self)
        if self.nGrowFootprint is None and self.nInitialRadii is None:
            raise ValueError("At least one of nGrowFootprint and nInitialRadii must be not None.")


def setupFitRegion(config, mask=None, detectionFootprint=None, modelEllipse=None, psfMoments=None):
    """!Return a Footprint containing the pixels that should be used in a model fit.

    The fit region is computed from the union of two regions:
     - the detection Footprint attached to the SourceRecord, grown by 'nGrowFootprint' pixels;
     - an initial model ellipse, scaled by 'nInitialRadii', then grown by the RMS size of the PSF.

    If either config parameter is None, that component of the fit region will be ignored entirely
    (the other component will be used directly, with no union).

    @param[in]  config    FitRegionConfig instance that controls how the region is determined.
    @param[in]  mask      afw.image.MaskU containing pixel mask information that could cause
                          pixels to be rejected from the fit (via config.badMaskPlanes).  If
                          None, config.badMaskPlanes is ignored.
    @param[in]  detectionFootprint    afw.detection. Footprint that will be grown by config.nGrowFootprint
                                      before being unioned into the fit region.  If None, the fit region
                                      will be set by the model ellipse region alone.
    @param[in]  modelEllipse    afw.geom.ellipses.Ellipse that will be scaled by config.nInitialRadii before
                                being unioned into the fit region.  If None, the fit region will be set by
                                the detection footprint region alone.
    @param[in]  psfMoments      afw.geom.ellipses.BaseCore instance that will be "convolved" (added in
                                quadrature) with the scaled model ellipse.  If None, the model ellipse
                                will be left unchanged after scaling.
    """
    if config.nGrowFootprint is not None and detectionFootprint is not None:
        if config.nGrowFootprint < 0:
            region1 = lsst.afw.detection.shrinkFootprint(detectionFootprint, -config.nGrowFootprint)
        elif config.nGrowFootprint > 0:
            region1 = lsst.afw.detection.growFootprint(detectionFootprint, config.nGrowFootprint)
        else:
            region1 = detectionFootprint
    else:
        region1 = None
    if config.nInitialRadii is not None and modelEllipse is not None:
        ellipse = lsst.afw.geom.ellipses.Ellipse(modelEllipse)
        ellipse.getCore().scale(config.nInitialRadii)
        if psfMoments is not None:
            ellipse = ellipse.convolve(lsst.afw.geom.ellipses.Ellipse(psfMoments))
        region2 = lsst.afw.detection.Footprint(ellipse, mask.getBBox(lsst.afw.image.PARENT))
    else:
        region2 = None
    if region1 and region2:
        region = lsst.afw.detection.mergeFootprints(region1, region2)
    elif region1:
        region = region1
    elif region2:
        region = region2
    else:
        raise ValueError("Neither detection Footprint nor modelEllipse specified for fit region.")
    region.clipTo(mask.getBBox(lsst.afw.image.PARENT))
    if mask is not None and config.badMaskPlanes:
        badPixelMask = 0x0
        for plane in config.badMaskPlanes:
            badPixelMask |= lsst.afw.image.MaskU.getPlaneBitMask(plane);
        region.intersectMask(mask, badPixelMask)
    return region

setupFitRegion.ConfigClass = FitRegionConfig
