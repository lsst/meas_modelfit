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
   python backgroundEstimationStageTest.py
or
   python
   >>> import backgroundEstimationStageTest
   >>> backgroundEstimationStageTest.run()
"""

import sys, os, math
from math import *

import unittest

import eups
import lsst.utils.tests as utilsTests
import lsst.pex.harness.Clipboard as pexClipboard
import lsst.pex.policy as pexPolicy
import lsst.meas.pipeline as measPipe
import lsst.afw.image as afwImage
import lsst.afw.display.ds9 as ds9

from lsst.pex.harness.simpleStageTester import SimpleStageTester

try:
    type(display)
except NameError:
    display = False

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class BackgroundEstimationStageTestCase(unittest.TestCase):
    """A test case for BackgroundEstimationStage.py"""

    def setUp(self):
        filename = os.path.join(eups.productDir("afwdata"), "CFHT", "D4", "cal-53535-i-797722_1")
        bbox = afwImage.BBox(afwImage.PointI(32,32), 512, 512)
        self.exposure = afwImage.ExposureF(filename, 0,bbox)        

    def tearDown(self):
        del self.exposure        

    def testSingleExposure(self):
        file = pexPolicy.DefaultPolicyFile("meas_pipeline", 
                                           "backgroundEstimation_policy.paf", "tests")
        policy = pexPolicy.Policy.createPolicy(file)

        stage = measPipe.BackgroundEstimationStage(policy)
        tester = SimpleStageTester(stage)

        clipboard = pexClipboard.Clipboard()         
        clipboard.put(policy.get("inputKeys.exposure"), self.exposure)

        if display:
            ds9.mtv(self.exposure, frame=0, title="Input")
        #
        # Do the work
        #
        outWorker = tester.runWorker(clipboard)

        outPolicy = policy.get("outputKeys")
        assert(outWorker.contains(outPolicy.get("backgroundSubtractedExposure")))
        assert(outWorker.contains(outPolicy.get("background")))

        if display:
            ds9.mtv(outWorker.get(outPolicy.get("backgroundSubtractedExposure")),
                    frame=1, title="Subtracted")

def suite():
    """Returns a suite containing all the test cases in this module."""

    utilsTests.init()

    suites = []

    if not eups.productDir("afwdata"):
        print >> sys.stderr, "afwdata is not setting up; skipping test"
    else:        
        suites += unittest.makeSuite(BackgroundEstimationStageTestCase)

    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)

