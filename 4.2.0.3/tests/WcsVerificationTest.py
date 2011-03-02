#!/usr/bin/env python

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

"""
Run with:
   python DetectTest.py

    This test isn't written yet, this is just boilerplate.
"""

import sys, os, math
from math import *

import pdb
import random

import unittest

import eups
import lsst.utils.tests as utilsTests
import lsst.pex.harness.Clipboard as pexClipboard
import lsst.pex.policy as pexPolicy
from lsst.pex.logging import Trace
import lsst.meas.pipeline as measPipe
import lsst.afw.image as afwImg
import lsst.afw.detection as afwDet


from lsst.pex.harness.simpleStageTester import SimpleStageTester

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def readSourceSet(fileName):
    fd = open(fileName, "r")

    sourceSet = afwDet.SourceSet()
    lineno = 0
    for line in fd.readlines():
        lineno += 1
        try:
            id, x, y, ra, dec, cts, flags = line.split()
        except Exception, e:
            print "Line %d: %s: %s" % (lineno, e, line),

        s = afwDet.Source()
        sourceSet.append(s)

        s.setId(int(id))
        s.setFlagForDetection(int(flags))
        s.setRa(float(ra))
        s.setXAstrom(float(x))
        s.setYAstrom(float(y))
        s.setDec(float(dec))
        s.setPsfFlux(float(cts))

    return sourceSet

class WcsVerificationStageTestCase(unittest.TestCase):
    """A test case for PhotoCalStage.py"""

    def setUp(self):
        random.seed(1369)
        
        ##Load sample input from disk
        path = os.path.join(eups.productDir("meas_pipeline"), "tests")
        srcSet = readSourceSet(os.path.join(path, "v695833-e0-c000.xy.txt"))
        
        #Make a copy, with different positions
        #The exact choice doesn't matter ,we just want to make sure the code returns an answer
        catSet = []
        for s in srcSet:
            s1 = afwDet.Source(s)
            
            s1.setXAstrom( s1.getXAstrom()+ random.uniform(-.1,.1))
            catSet.append(s1)
            
        #Make a SourceMatch object
        maxDist = 1/3600. #matches must be this close together
        srcMatchSet = afwDet.matchXy(catSet, srcSet, maxDist)
        
        #Put them on the clipboard
        filename = pexPolicy.DefaultPolicyFile("meas_pipeline", 
                      "WcsVerificationStageDictionary.paf", "policy")
        self.policy = pexPolicy.Policy.createPolicy(filename)

        self.clipboard = pexClipboard.Clipboard()         
        self.clipboard.put(self.policy.get("sourceMatchSetKey"), srcMatchSet)
        self.clipboard.put(self.policy.get("inputExposureKey"),
                afwImg.ExposureF(4000, 4000))

    def tearDown(self):
        del self.clipboard
        del self.policy
        
    def testSingleExposure(self):
        #Run the stage
        stage = measPipe.WcsVerificationStage(self.policy)
        tester = SimpleStageTester(stage)
        outWorker = tester.runWorker(self.clipboard)

        #Check for output parameters
        image = self.clipboard.get(self.policy.get("inputExposureKey"))
        metadata = image.getMetadata()
        self.assert_(metadata.exists("minObjectsPerCell"))
        self.assert_(metadata.exists("maxObjectsPerCell"))
        self.assert_(metadata.exists("meanObjectsPerCell"))
        self.assert_(metadata.exists("stdObjectsPerCell"))

    
def suite():
    """Returns a suite containing all the test cases in this module."""

    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(WcsVerificationStageTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)

