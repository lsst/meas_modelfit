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
from lsst.pex.harness.simpleStageTester import SimpleStageTester
import lsst.pex.policy as pexPolicy
from lsst.pex.logging import Log, Debug, LogRec, Prop
from lsst.pex.exceptions import LsstCppException
import lsst.afw.image as afwImg

import lsst.meas.astrom as measAstrom
import lsst.meas.astrom.net as astromNet
import lsst.meas.astrom.sip as sip
import lsst.meas.photocal as photocal
import lsst.meas.astrom.sip.cleanBadPoints as cleanBadPoints

import pdb
    

class PhotoCalStageParallel(harnessStage.ParallelProcessing):
    """Calculate the magnitude zero point for a SourceSet for an image that
    has been matched to a corresponding SourceSet for a catalogue
    """
    
    def setup(self):
        policyFile=pexPolicy.DefaultPolicyFile("meas_pipeline",   # package name
                                  "PhotoCalStageDictionary.paf", # default. policy
                                  "policy" # dir containing policies
                                  )
        defaultPolicy = pexPolicy.Policy.createPolicy(policyFile, policyFile.getRepositoryPath())
        
        if self.policy is None:
            self.policy = defaultPolicy
        else:
            self.policy.mergeDefaults(defaultPolicy)
               
        #Setup the log
        self.log = Debug(self.log, "PhotoCalStageParallel")
        self.log.setThreshold(Log.DEBUG)
        self.log.log(Log.INFO, "Finished setup of PhotoCalStageParallel")
        

    def process(self, clipboard):
        self.log.log(Log.INFO, "Determining Photometric Zero Point")

        #Check inputs
        if clipboard is None:
            raise RuntimeError("Clipboard is empty")

        expKey = self.policy.get('inputExposureKey')
        if not clipboard.contains(expKey):
            raise RuntimeError("No exposure on clipboard")
        exp = clipboard.get(expKey)

        srcMatchSetKey = self.policy.get("sourceMatchSetKey")
        if not clipboard.contains(srcMatchSetKey):
            raise RuntimeError("No input SourceMatch set on clipboard")
        srcMatchSet = clipboard.get(srcMatchSetKey)            
        
       
        
        #Do the work
        try:
            magObj = photocal.calcPhotoCal(srcMatchSet, log=self.log)
        except ValueError, e:
            msg = "Failed to calculate photometric zeropoint: %s" %(e)
            self.log.log(Log.FATAL, msg)
            magObj = None

        if magObj is not None:
            exp.getCalib().setFluxMag0(magObj.getFlux(0))
            self.log.log(Log.INFO, "Flux of magnitude 0: %g" %
                    (magObj.getFlux(0),))

        #Save results to clipboard
        outputValueKey = self.policy.get("outputValueKey")
        clipboard.put(outputValueKey, magObj)



class PhotoCalStage(harnessStage.Stage):
    """A wrapper stage that supplies the names of the classes that do the work
       Different classes are provided for serial and parallel processing
    """
    parallelClass = PhotoCalStageParallel




