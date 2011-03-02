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

class SourceDetectionStageParallel(harnessStage.ParallelProcessing):
    """
    Description:
        This stage wraps the detection of sources on an exposure provided
        on the clipboard. If more than one exposure are provided, the exposures
        will be added first to create a stacked exposure.
        

    Policy Dictionaty:
    lsst/meas/pipeline/SourceDetectionStageDictionary.paf

    Clipboard Input:
    - Calibrated science Exposure(s) (including background)
    - a PSF may be specified by policy attribute inputPsf. Alternatively, the
      stage's policy may request that a psf be constructed, by providing the
      psfPolicy attribute.

    ClipboardOutput:
    - background subtracted Exposure used in the detection. Key specified
        by policy attribute 'backgroundSubtractedExposure'
    - the measured background object itself. Key specified by policy 
        attribute 'background'        
    - PSF: the psf used to smooth the exposure before detection 
        Key specified by policy attribute 'psfKey'
    - PositiveFootprintSet (DetectionSet): if thresholdPolarity policy 
        is "positive" or "both". Key specified by policy attribute
        'positiveDetectionKey'
    - NegativeFootprintSet (DetectionSet): if threholdPolarity policy 
        is "negative" or "both". Key specified by policy attribute
        'negativeDetectionKey'
    """
    def setup(self):
        self.log = Log(self.log, "SourceDetectionStage - parallel")

        policyFile = pexPolicy.DefaultPolicyFile("meas_pipeline", 
            "SourceDetectionStageDictionary.paf", "policy")
        defPolicy = pexPolicy.Policy.createPolicy(policyFile, policyFile.getRepositoryPath(), True)

        if self.policy is None:
            self.policy = pexPolicy.Policy()
        self.policy.mergeDefaults(defPolicy.getDictionary())

    def process(self, clipboard):
        """
        Detect sources in the worker process
        """
        self.log.log(Log.INFO, "Detecting Sources in process")
        
        #grab exposure(s) from clipboard
        exposureList = []        
        if self.policy.isArray("inputKeys.exposure"):
            #if more than one exposure is supplied, first
            exposureKeyList = self.policy.getArray("inputKeys.exposure")
            for key in exposureKeyList:
                exposureList.append(clipboard.get(key))                        
    
            #stack the exposures
            exposure = sourceDetection.addExposures(exposureList)
        else:                    
            exposure = clipboard.get(self.policy.get("inputKeys.exposure"))
            exposureList.append(exposure)
       
        if exposure is None:
            raise RuntimeError, "Input Exposure Missing"
            
        #subtract the background
        if self.policy.exists("backgroundPolicy") and self.policy.get("backgroundPolicy.algorithm") != "NONE":
            background, backgroundSubtracted = sourceDetection.estimateBackground(
                exposure, self.policy.get("backgroundPolicy"), True)
        else:
            backgroundSubtracted, background = exposure, None

        #get or make a smoothing psf according to the policy
        if self.policy.exists("inputKeys.psf"):
            psf = clipboard.get(self.policy.get("inputKeys.psf"))
        elif self.policy.exists("psfPolicy"):
            psf = sourceDetection.makePsf(self.policy.get("psfPolicy"))
        else:
            psf = None
            
        #perform detection
        dsPositive, dsNegative = sourceDetection.detectSources(
            backgroundSubtracted, psf, self.policy.get("detectionPolicy"))        

        #
        # Copy detection mask bits to the input exposure(s)
        #
        detectionMask = backgroundSubtracted.getMaskedImage().getMask()
        for e in exposureList:
            msk = e.getMaskedImage().getMask()
            msk |= detectionMask
            del msk

        del detectionMask


        #output products
        if dsPositive:
            clipboard.put(
                self.policy.get("outputKeys.positiveDetection"), dsPositive)
        if dsNegative:
            clipboard.put(
                self.policy.get("outputKeys.negativeDetection"), dsNegative)
        if background:
            clipboard.put(self.policy.get("outputKeys.background"), background)
        if psf:
            clipboard.put(self.policy.get("outputKeys.psf"), psf)

        clipboard.put(self.policy.get("outputKeys.backgroundSubtractedExposure"),
                      backgroundSubtracted)
        
class SourceDetectionStage(harnessStage.Stage):
    parallelClass = SourceDetectionStageParallel

