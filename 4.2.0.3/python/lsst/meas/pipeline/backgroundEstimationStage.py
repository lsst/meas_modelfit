#! python

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

import lsst.pex.harness.stage as harnessStage

from lsst.pex.logging import Log

import lsst.pex.policy as pexPolicy
import lsst.afw.detection as afwDet
import lsst.afw.image as afwImg
import lsst.afw.math as afwMath
import lsst.pex.exceptions as pexExcept
import lsst.meas.algorithms as measAlg

import lsst.meas.utils.sourceDetection as sourceDetection

import lsst.afw.display.ds9 as ds9
import lsst.afw.display.utils as displayUtils

try:
    type(display)
except NameError:
    display = False

class BackgroundEstimationStageParallel(harnessStage.ParallelProcessing):
    """
    Description:
        This stage wraps estimating and possibly subtracting the background from an exposure
        on the clipboard.        

    Policy Dictionary:
    lsst/meas/pipeline/policy/BackgroundEstimationStageDictionary.paf

    Clipboard Input:
    - Calibrated science Exposure(s) (including background)

    ClipboardOutput:
    - background subtracted Exposure used in the detection. Key specified
        by policy attribute 'backgroundSubtractedExposure'
    - the measured background object itself. Key specified by policy 
        attribute 'background'        
    """
    def setup(self):
        self.log = Log(self.log, "BackgroundEstimationStage - parallel")

        policyFile = pexPolicy.DefaultPolicyFile("meas_pipeline", 
                                                 "BackgroundEstimationStageDictionary.paf", "policy")
        defPolicy = pexPolicy.Policy.createPolicy(policyFile, policyFile.getRepositoryPath(), True)

        if self.policy is None:
            self.policy = pexPolicy.Policy()
        self.policy.mergeDefaults(defPolicy.getDictionary())

    def process(self, clipboard):
        """
        Detect sources in the worker process
        """
        self.log.log(Log.INFO, "Subtracting background in process")
        
        #grab exposure from clipboard
        exposure = clipboard.get(self.policy.get("inputKeys.exposure"))
            
        #estimate and maybe subtract the background
        background, backgroundSubtractedExposure = sourceDetection.estimateBackground(
            exposure,
            self.policy.get("parameters.backgroundPolicy"),
            self.policy.get("parameters.subtractBackground"))

        #output products
        clipboard.put(self.policy.get("outputKeys.background"), background)
        if backgroundSubtractedExposure:
            clipboard.put(self.policy.get("outputKeys.backgroundSubtractedExposure"),
                          backgroundSubtractedExposure)
        
        
class BackgroundEstimationStage(harnessStage.Stage):
    parallelClass = BackgroundEstimationStageParallel

