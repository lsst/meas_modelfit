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

import lsst.pex.harness as pexHarness
import lsst.pex.harness.stage as harnessStage
from lsst.pex.harness.simpleStageTester import SimpleStageTester
import lsst.pex.policy as pexPolicy
from lsst.pex.logging import Log, Debug, LogRec, Prop
from lsst.pex.exceptions import LsstCppException
import lsst.afw.image as afwImg
import lsst.afw.detection as afwDet

import lsst.meas.astrom as measAstrom
import lsst.meas.astrom.net as astromNet
import lsst.meas.astrom.sip as sip
import lsst.meas.astrom.sip.cleanBadPoints as cleanBadPoints

class WcsDeterminationStageParallel(harnessStage.ParallelProcessing):
    """Validate the Wcs for an image using the astrometry.net package and calculate distortion
    coefficients
    
    Given a initial Wcs, and a list of sources (with pixel positions for each) in an image,
    pass these to the astrometry_net package to verify the result. Then calculate
    the distortion in the image and add that to the Wcs as SIP polynomials
    
    Clipboard Input:
    - an Exposure containing a Wcs
    - a SourceSet
    
    Clipboard Output
    - A wcs
    - A vector of SourceMatch objects
    """
    
    def setup(self):
        #I don't have the default policy in the correct place yet
        policyFile=pexPolicy.DefaultPolicyFile("meas_pipeline",   # package name
                                  "WcsDeterminationStageDictionary.paf", # default. policy
                                  "policy" # dir containing policies
                                  )
        defaultPolicy = pexPolicy.Policy.createPolicy(policyFile, policyFile.getRepositoryPath())
        
        #The stage can be called with an optional local policy file, which overrides the defaults
        #merge defaults
        policyFile = pexPolicy.DefaultPolicyFile("meas_pipeline",
            "WcsDeterminationStageDictionary.paf", "policy")
        defPolicy = pexPolicy.Policy.createPolicy(policyFile, policyFile.getRepositoryPath(), True)

        if self.policy is None:
            self.policy = defaultPolicy
        else:
            self.policy.mergeDefaults(defaultPolicy)


        #Setup the astrometry solver
        path=os.path.join(os.environ['ASTROMETRY_NET_DATA_DIR'], "metadata.paf")
        self.solver = astromNet.GlobalAstrometrySolution(path)
        # self.solver.allowDistortion(self.policy.get('allowDistortion'))
        self.solver.setMatchThreshold(self.policy.get('matchThreshold'))
               
        #Setup the log
        self.log = Debug(self.log, "WcsDeterminationStageParallel")
        self.log.setThreshold(Log.DEBUG)
        self.log.log(Log.INFO, "Finished setup of WcsDeterminationStageParallel")
        

    def process(self, clipboard):
        self.log.log(Log.INFO, "Determining Wcs")

        #Check inputs
        if clipboard is None:
            raise RuntimeError("Clipboard is empty")

        expKey = self.policy.get('inputExposureKey')
        if not clipboard.contains(expKey):
            raise RuntimeError("No exposure on clipboard")
        exp = clipboard.get(expKey)

        srcSetKey=self.policy.get('inputSourceSetKey')
        
        if not clipboard.contains(srcSetKey):
            raise RuntimeError("No wcsSourceSet on clipboard")
        srcSet = clipboard.get(srcSetKey)
        
        # Determine list of matching sources, and WCS
        astrom  = measAstrom.determineWcs(self.policy, exp, 
                                          srcSet, solver=self.solver, log=self.log)
        matchList = astrom.getMatches()
        wcs = astrom.getWcs()
        matchListMeta = astrom.getMatchMetadata()

        #Save results to clipboard
        smv = afwDet.SourceMatchVector()
        for m in matchList:
            smv.push_back(m)
        psmv = afwDet.PersistableSourceMatchVector(smv, matchListMeta)

        clipboard.put(self.policy.get('outputMatchListKey'), matchList)
        clipboard.put(self.policy.get('outputMatchListMetaKey'), matchListMeta)
        clipboard.put(self.policy.get('outputMatchListKey') + '_persistable', psmv)
        clipboard.put(self.policy.get('outputWcsKey'), wcs)


class WcsDeterminationStage(harnessStage.Stage):
    """A wrapper stage that supplies the names of the classes that do the work
       Different classes are provided for serial and parallel processing
    """
    parallelClass = WcsDeterminationStageParallel




