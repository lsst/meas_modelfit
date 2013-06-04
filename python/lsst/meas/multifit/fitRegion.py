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
import lsst.afw.detection

__all__ = ("FitRegionConfig", "setupFitRegion")

class FitRegionConfig(lsst.pex.config.Config):
    """Config class for setupFitRegion() function
    """
    nGrow = lsst.pex.config.Field(
        dtype=int,
        default=5,
        doc="How many pixels to grow detection footprints by"
    )
    growIsotropic = lsst.pex.config.Field(
        dtype=bool,
        default=True,
        doc="Whether to grow footprints isotropically"
    )

def setupFitRegion(config, exposure, source):
    """Given a SourceRecord (with Footprint) and the Exposure it was detected on,
    return a new Footprint containing the pixels that should be used in a model
    fit of the given source.

    In the future, this might do something fancier than just grow the footprint, such as masking
    out bad pixels or ORing the detection footprint with an ellipse derived from the source, but
    we don't anticipate needing any of that for S13.
    """
    return lsst.afw.detection.growFootprint(source.getFootprint(), config.nGrow, config.growIsotropic)
setupFitRegion.ConfigClass = FitRegionConfig
