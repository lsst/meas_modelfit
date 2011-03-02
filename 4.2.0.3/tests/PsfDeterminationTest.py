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
   python PsfDeterminationTest.py
"""

import sys, os, math
from math import *

import pdb
import unittest

import eups
import lsst.utils.tests as utilsTests
import lsst.pex.harness.Clipboard as pexClipboard
import lsst.pex.policy as pexPolicy
import lsst.meas.pipeline as measPipe
import lsst.afw.detection as afwDet
import lsst.afw.image as afwImage
from lsst.pex.harness.simpleStageTester import SimpleStageTester

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class PsfDeterminationStageTestCase(unittest.TestCase):
    """A test case for SourceMeasurementStage.py"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def clipboardIoTest(self):
        file = pexPolicy.DefaultPolicyFile("meas_pipeline", 
                "tests/sourceDetection0_policy.paf")
        detectPolicy = pexPolicy.Policy.createPolicy(file)
        
        file = pexPolicy.DefaultPolicyFile("meas_pipeline", 
                "tests/sourceMeasurement0_policy.paf")
        measurePolicy = pexPolicy.Policy.createPolicy(file)

        file = pexPolicy.DefaultPolicyFile("meas_pipeline", 
                "tests/psfDetermination0_policy.paf")
        psfPolicy = pexPolicy.Policy.createPolicy(file)
       

        tester = SimpleStageTester()
        tester.addStage(measPipe.SourceDetectionStage(detectPolicy))
        tester.addStage(measPipe.SourceMeasurementStage(measurePolicy))
        tester.addStage(measPipe.PsfDeterminationStage(psfPolicy))


        clipboard = pexClipboard.Clipboard() 
        filename = os.path.join(eups.productDir("afwdata"),
                                "CFHT", "D4", 
                                "cal-53535-i-797722_1")
        # test only a portion of the exposure to speed up testing
        bbox = afwImage.BBox(afwImage.PointI(32, 32), 512, 512)        
        testExp = afwImage.ExposureF(filename, 0, bbox)

        # test full image
        # testExp = afImage.ExposureF(filename)

        clipboard = pexClipboard.Clipboard() 
        clipboard.put(detectPolicy.get("inputKeys.exposure"), testExp)
        
        
        outWorker = tester.runWorker(clipboard)
  
        assert(outWorker.contains(psfPolicy.get("data.outputPsfKey")))
        assert(outWorker.contains(psfPolicy.get("data.outputCellSetKey")))

        del outWorker
        del testExp
        del clipboard
        del tester

        print >> sys.stderr, "at end of test"

def suite():
    """Returns a suite containing all the test cases in this module."""

    utilsTests.init()

    suites = []
    if not eups.productDir("afwdata"):
        print >> sys.stderr, "afwdata is not setting up; skipping test"
    else:
        suites += unittest.makeSuite(PsfDeterminationStageTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)

