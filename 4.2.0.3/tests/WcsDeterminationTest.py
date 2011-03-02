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
"""

import sys, os, math
from math import *

import pdb
import unittest

import eups
import lsst.utils.tests as utilsTests
import lsst.pex.harness.Clipboard as pexClipboard
import lsst.pex.policy as pexPolicy
from lsst.pex.logging import Trace
import lsst.meas.pipeline as measPipe
import lsst.afw.detection as afwDet
import lsst.afw.image as afwImg


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

class WcsDeterminationStageTestCase(unittest.TestCase):
    """A test case for WcsDeterminationStage.py"""

    def setUp(self):
        #Load sample input from disk
        path = os.path.join(eups.productDir("meas_pipeline"), "tests")
        exp = afwImg.ExposureF(os.path.join(path, "v695833-e0-c000-a00.sci"))
        srcSet = readSourceSet(os.path.join(path, "v695833-e0-c000.xy.txt"))

        #Put them on the clipboard
        fileName = pexPolicy.DefaultPolicyFile("meas_pipeline", 
                "WcsDeterminationStageDictionary.paf", "policy")
        self.policy = pexPolicy.Policy.createPolicy(fileName)

        self.clipboard = pexClipboard.Clipboard()         
        self.clipboard.put(self.policy.get("inputExposureKey"), exp)
        self.clipboard.put(self.policy.get("inputSourceSetKey"), srcSet)

        
    def tearDown(self):
        del self.policy
        del self.clipboard

    def testSingleExposure(self):
        #Run the stage
        stage = measPipe.WcsDeterminationStage(self.policy)
        tester = SimpleStageTester(stage)
        outWorker = tester.runWorker(self.clipboard)

        #Check for output parameters
        wcsKey = self.policy.get("outputWcsKey")
        assert(outWorker.contains(wcsKey))

        matchListKey = self.policy.get("outputMatchListKey")
        assert(outWorker.contains(matchListKey))
        assert(len(outWorker['matchListKey']) > 0)


        
def suite():
    """Returns a suite containing all the test cases in this module."""

    utilsTests.init()

    suites = []

    if not eups.productDir("astrometry_net_data"):
        print >> sys.stderr, "Unable to test WcsDeterminationStage as astrometry_net_data is not setup"
    else:
        suites += unittest.makeSuite(WcsDeterminationStageTestCase)

    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)

