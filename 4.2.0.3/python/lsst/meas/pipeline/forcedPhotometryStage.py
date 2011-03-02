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

from lsst.pex.logging import Log
import lsst.pex.harness.stage as harnessStage
import lsst.pex.policy as pexPolicy

__all__ = ["ForcedPhotometryStage", "ForcedPhotometryStageParallel"]

class ForcedPhotometryStageParallel(harnessStage.ParallelProcessing):
    """
    Given model(s) for an object, measure the source on the exposures
    provided in the ExposureStack

    Clipboard input:
    - Model(s): point source and/or small galaxy model    
    - ExposureStack
    """

    def setup(self):
        self.log = Log(self.log, "ForcedPhotometryStage - parallel")

        policyFile = pexPolicy.DefaultPolicyFile("meas_pipeline", 
            "ForcedPhotometryStageDictionary.paf", "policy")
        defPolicy = pexPolicy.Policy.createPolicy(policyFile, policyFile.getRepositoryPath(), True)

        if self.policy is None:
            self.policy = pexPolicy.Policy()
        self.policy.mergeDefaults(defPolicy.getDictionary())

    def process(self, clipboard):
        pass

class ForcedPhotometryStage(harnessStage.Stage):
    parallelClass = ForcedPhotometryStageParallel
