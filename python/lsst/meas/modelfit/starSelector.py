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
import lsst.meas.algorithms
from lsst.meas.algorithms.starSelectorRegistry import starSelectorRegistry

class S13StarSelectorConfig(lsst.pex.config.Config):
    fluxMin = lsst.pex.config.Field(
        doc = "specify the minimum apFlux for good PsfCandidates",
        dtype = float,
        default = 50000,
        check = lambda x: x >= 0.0,
    )
    kernelSize = lsst.pex.config.Field(
        doc = "size of the Psf kernel to create",
        dtype = int,
        default = 21,
    )
    borderWidth = lsst.pex.config.Field(
        doc = "number of pixels to ignore around the edge of PSF candidate postage stamps",
        dtype = int,
        default = 0,
    )

class S13StarSelector(object):
    ConfigClass = S13StarSelectorConfig

    def __init__(self, config):
        self.config = config

    def selectStars(self, exposure, catalog, matches=None):
        psfCandidateList = []
        for source in catalog:
            if source.getApFluxFlag():
                continue
            if source.getApFlux() < self.config.fluxMin:
                continue
            psfCandidate = lsst.meas.algorithms.makePsfCandidate(source, exposure)
            if psfCandidate.getWidth() == 0:
                psfCandidate.setBorderWidth(self.config.borderWidth)
                psfCandidate.setWidth(self.config.kernelSize + 2*self.config.borderWidth)
                psfCandidate.setHeight(self.config.kernelSize + 2*self.config.borderWidth)
            psfCandidateList.append(psfCandidate)
        return psfCandidateList

starSelectorRegistry.register("s13", S13StarSelector)
