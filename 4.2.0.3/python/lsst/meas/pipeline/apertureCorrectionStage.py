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
from math import *

from lsst.pex.logging import Log
import lsst.pex.harness.stage as harnessStage
import lsst.pex.policy as pexPolicy
import lsst.meas.algorithms as measAlg
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.meas.algorithms.ApertureCorrection as apertureCorrection
import lsst.sdqa as sdqa

class ApertureCorrectionStageParallel(harnessStage.ParallelProcessing):
    """
    Given an exposure and a set of sources measured on that exposure,
    determine the aperture correction for that exposure.

    This stage works on lists of (exposure, sourceSet) pairs.

    Their location on the clipboard is specified via policy.
    see lsst/meas/pipeline/pipeline/ApertureCorrectionStageDictionary.paf
    for details on configuring valid stage policies
    """
    def setup(self):
        self.log = Log(self.log, "ApertureCorrectionStage - parallel")

        # aperture correction policy
        apCorrPolicyFile = pexPolicy.DefaultPolicyFile("meas_pipeline", 
                                                       "ApertureCorrectionStageDictionary.paf", "policy")
        defPolicy = pexPolicy.Policy.createPolicy(apCorrPolicyFile,
                                                  apCorrPolicyFile.getRepositoryPath(), True)

        if self.policy is None:
            self.policy = pexPolicy.Policy()
        self.policy.mergeDefaults(defPolicy.getDictionary())

        self.ApCorrPolicy = self.policy.get("parameters.ApertureCorrectionPolicy")

    def process(self, clipboard):
        self.log.log(Log.INFO, "Estimating Aperture Correction is in process")

        
        #grab exposure from clipboard
        exposure = clipboard.get(self.policy.get("inputKeys.exposure"))       
        cellSet = clipboard.get(self.policy.get("inputKeys.cellSet"))
        
        sdqaRatings = sdqa.SdqaRatingSet()
        apCorrCtrl = apertureCorrection.ApertureCorrectionControl(self.ApCorrPolicy)
        apCorr = apertureCorrection.ApertureCorrection(exposure, cellSet, sdqaRatings,
                                                       apCorrCtrl, log=self.log)
        

        clipboard.put(self.policy.get("outputKeys.apCorr"), apCorr)
        clipboard.put(self.policy.get("outputKeys.sdqa"), sdqa)

        
class ApertureCorrectionStage(harnessStage.Stage):
    parallelClass = ApertureCorrectionStageParallel

