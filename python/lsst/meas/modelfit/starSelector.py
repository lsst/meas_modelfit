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
from lsst.afw.table import SourceCatalog
from lsst.pipe.base import Struct
from lsst.meas.algorithms import StarSelectorTask, starSelectorRegistry

__all__ = ["S13StarSelectorConfig", "S13StarSelectorTask"]

class S13StarSelectorConfig(StarSelectorTask.ConfigClass):
    fluxMin = lsst.pex.config.Field(
        doc = "specify the minimum apFlux for good PsfCandidates",
        dtype = float,
        default = 50000,
        check = lambda x: x >= 0.0,
    )

class S13StarSelectorTask(StarSelectorTask):
    ConfigClass = S13StarSelectorConfig
    usesMatches = False # selectStars does not use its matches argument

    def selectStars(self, exposure, sourceCat, matches=None):
        starCat = SourceCatalog(sourceCat.schema)
        for source in sourceCat:
            if source.getApFluxFlag():
                continue
            if source.getApFlux() < self.config.fluxMin:
                continue
            starCat.append(source)
        return Struct(
            starCat = starCat,
        )

starSelectorRegistry.register("s13", S13StarSelectorTask)
