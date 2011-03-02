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


import lsst.pex.harness.stage as harnessStage
from lsst.pex.logging import Log
import lsst.daf.base as dafBase
from lsst.daf.base import *
import lsst.pex.policy as pexPolicy

class TransformDetectionStageParallel(harnessStage.ParallelProcessing):
    """
    This stage takes a single Model. 
    Computes an ra/dec bounding box for that model.
    Uses the ImageAccessApi to pull down a list of images
      that overlap that bounding box
    For each image in the list, 
        Uses the ImageAccessApi to pull down the exposure's metadata
        computes the model's projection's bbox on that exposure
          using the metadata
        adds the image filename and bbox to a PropertrySet
    puts the property set on the clipboard for future stages to use

    INPUT:
    - a Model
    - WCS list
    - Psf list
    OUTPUT:
    - BBox list
    """
    
    def setup(self):
        self.log = Log(self.log, "TransformDetectionStage - parallel")
 
        policyFile = pexPolicy.DefaultPolicyFile("meas_pipeline", 
            "TransformDetectionStageDictionary.paf", "policy")
        defPolicy = pexPolicy.Policy.createPolicy(policyFile, policyFile.getRepositoryPath(), True)

        if self.policy is None:
            self.policy = pexPolicy.Policy()
        self.policy.mergeDefaults(defPolicy.getDictionary())        

    def process(self, clipboard):
        model = clipboard.get(self.policy.get("inputKeys.initialSGModel"))
        quadsphere = quadsphere.makeQuadsphere(self.policy.get("parameters.quadspherePolicy"))
        pixIds = measUtils.multifit.computeSkypixIdList(model, quadsphere)
        clipboard.put(self.policy.get("outputKeys.skypixId"), pixIds)


