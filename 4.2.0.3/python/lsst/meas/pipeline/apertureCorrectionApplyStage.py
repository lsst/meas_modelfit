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
import lsst.afw.detection as afwDet
import lsst.meas.algorithms.ApertureCorrection as apertureCorrection
import lsst.sdqa as sdqa

class ApertureCorrectionApplyStageParallel(harnessStage.ParallelProcessing):
    """
    Given a set of sources measured, and an apertureCorrection object,
    apply the aperture correction to the sourceSet.

    Their location on the clipboard is specified via policy.
    see lsst/meas/pipeline/pipeline/ApertureCorrectionApplyStageDictionary.paf
    for details on configuring valid stage policies
    """
    def setup(self):
        self.log = Log(self.log, "ApertureCorrectionApplyStage - parallel")

        # aperture correction policy
        apCorrPolicyFile = pexPolicy.DefaultPolicyFile("meas_pipeline", 
                                                       "ApertureCorrectionApplyStageDictionary.paf", "policy")
        defPolicy = pexPolicy.Policy.createPolicy(apCorrPolicyFile,
                                                  apCorrPolicyFile.getRepositoryPath(), True)
        
        if self.policy is None:
            self.policy = pexPolicy.Policy()
        self.policy.mergeDefaults(defPolicy.getDictionary())

    def process(self, clipboard):
        self.log.log(Log.INFO, "Aperture Correction Apply is in process")

        
        #grab sourceSet and apertureCorrection from clipboard
        # correct psf flux in situ
        sourceSet = clipboard.get(self.policy.get("inputKeys.sourceSet"))
        apCorr    = clipboard.get(self.policy.get("inputKeys.apCorr"))

        for s in sourceSet:
            
            # apply the correction, and propegate the error
            ac, acErr = apCorr.computeAt(s.getXAstrom(), s.getYAstrom())
            s.setPsfFlux(s.getPsfFlux()*ac)
            varFlux = ac**2 * s.getPsfFluxErr()**2
            varApCorr = s.getPsfFlux()**2 * acErr**2
            fluxErrPropegated = math.sqrt(varFlux + varApCorr) # assume covariance = 0
            s.setPsfFluxErr(fluxErrPropegated)

        
class ApertureCorrectionApplyStage(harnessStage.Stage):
    parallelClass = ApertureCorrectionApplyStageParallel

