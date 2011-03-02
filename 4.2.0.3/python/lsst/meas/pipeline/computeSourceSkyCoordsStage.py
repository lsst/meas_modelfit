# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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

import math

import lsst.pex.harness.stage as harnessStage

from lsst.pex.logging import Log

import lsst.pex.policy as pexPolicy
import lsst.afw.detection as afwDet
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.pex.exceptions as pexExcept
import lsst.meas.utils.sourceMeasurement as sourceMeasurement

class ComputeSourceSkyCoordsStageParallel(harnessStage.ParallelProcessing):
    """
    Description:
       Stage that converts pixel coordinates (X,Y) to sky coordinates
       (ra,dec) (in radians) for all sources in a SourceSet. Note
       that sources are updated in place and that the clipboard structure
       is not changed in any way.

    Policy Dictionary:
        lsst/meas/pipeline/SourceXYToRaDecStageDictionary.paf
    """
    def setup(self):
        self.log = Log(self.log, "ComputeSourceSkyCoordsStage - parallel")
        policyFile = pexPolicy.DefaultPolicyFile(
            "meas_pipeline", "ComputeSourceSkyCoordsStageDictionary.paf", "policy")
        defPolicy = pexPolicy.Policy.createPolicy(
            policyFile, policyFile.getRepositoryPath(), True)
        if self.policy is None:
            self.policy = pexPolicy.Policy()
        self.policy.mergeDefaults(defPolicy.getDictionary())

    def process(self, clipboard):
        wcsKey = self.policy.getString("inputKeys.wcs")
        sourcesKey = self.policy.getString("inputKeys.sources")
        exposureKey = self.policy.getString("inputKeys.exposure")
        if clipboard.contains(wcsKey):
            wcs = clipboard.get(wcsKey)
        else:
            exposure = clipboard.get(exposureKey)
            wcs = exposure.getWcs()
        sourceSet = clipboard.get(sourcesKey)
        sourceMeasurement.computeSkyCoords(wcs, sourceSet)
         

class ComputeSourceSkyCoordsStage(harnessStage.Stage):
    parallelClass = ComputeSourceSkyCoordsStageParallel
