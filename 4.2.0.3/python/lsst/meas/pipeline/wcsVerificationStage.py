#! /usr/bin/env python

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


import os
import sys
import math

import eups
import lsst.pex.harness as pexHarness
import lsst.pex.harness.stage as harnessStage
import lsst.pex.policy as pexPolicy
from lsst.pex.logging import Log, Debug, LogRec, Prop
from lsst.pex.exceptions import LsstCppException

import lsst.meas.astrom.sip as sip
import lsst.meas.astrom.verifyWcs as verifyWcs

class WcsVerificationParallel(harnessStage.ParallelProcessing):
    """Compute some statistics that indicate if we did a good job computing the Wcs for an image.
    """
    
    def setup(self):
        policyFile=pexPolicy.DefaultPolicyFile("meas_pipeline",   # package name
                                  "WcsVerificationStageDictionary.paf", # default. policy
                                  "policy" # dir containing policies
                                  )
        defaultPolicy = pexPolicy.Policy.createPolicy(policyFile, policyFile.getRepositoryPath())
        
        if self.policy is None:
            self.policy = defaultPolicy
        else:
            self.policy.mergeDefaults(defaultPolicy)
               
        #Setup the log
        self.log = Debug(self.log, "WcsVerificationStageParallel")
        self.log.setThreshold(Log.DEBUG)
        self.log.log(Log.INFO, "Finished setup of WcsVerificationStageParallel")
        

    def process(self, clipboard):
        self.log.log(Log.INFO, "Calculating statistics on wcs solution")

        #Get input
        if clipboard is None:
            raise RuntimeError("Clipboard is empty")

        srcMatchSetKey = self.policy.get("sourceMatchSetKey")
        if not clipboard.contains(srcMatchSetKey):
            raise RuntimeError("No input SourceMatch set on clipboard")
        srcMatchSet = clipboard.get(srcMatchSetKey)

        inputExposureKey = self.policy.get("inputExposureKey")
        if clipboard.contains(inputExposureKey):
            exposure = clipboard.get(inputExposureKey)
        else:
            exposure = None
        
        #Do the work
        outputDict = {}
        outputDict.update(sip.sourceMatchStatistics(srcMatchSet))
        outputDict.update(verifyWcs.checkMatches(srcMatchSet, exposure, self.log))

        self.log.log(self.log.DEBUG, "cells nobj min = %(minObjectsPerCell)s max = %(maxObjectsPerCell)s mean = %(meanObjectsPerCell)s std = %(stdObjectsPerCell)s" % outputDict)
        #
        # Set the metadata
        #
        for k, v in outputDict.items():
            exposure.getMetadata().set(k, v)

class WcsVerificationStage(harnessStage.Stage):
    """A wrapper stage that supplies the names of the classes that do the work
       Different classes are provided for serial and parallel processing
    """
    parallelClass = WcsVerificationParallel
