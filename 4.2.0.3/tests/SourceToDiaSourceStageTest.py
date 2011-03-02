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
   python SourceToDiaSourceStageTest.py
"""

import sys, os, math
from math import *

import pdb
import unittest
import random
import time

import lsst.utils.tests as utilsTests
import lsst.daf.base as dafBase
import lsst.pex.harness.Queue as pexQueue
import lsst.pex.harness.Clipboard as pexClipboard
import lsst.pex.policy as pexPolicy
import lsst.meas.pipeline as measPipe
import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
from lsst.pex.harness.simpleStageTester import SimpleStageTester

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

class SourceToDiaSourceStageTestCase(unittest.TestCase):
    """A test case for SourceDetectionStage.py"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testValidClipboard(self):
        dataPolicy1 = pexPolicy.Policy()
        dataPolicy1.add("inputKey", "sourceSet0")
        dataPolicy1.add("outputKey", "diaSourceSet0")
        dataPolicy2 = pexPolicy.Policy()
        dataPolicy2.add("inputKey", "sourceSet1")
        dataPolicy2.add("outputKey", "diaSourceSet1")
        stagePolicy = pexPolicy.Policy()
        stagePolicy.add("data", dataPolicy1)
        stagePolicy.add("data", dataPolicy2)
        stagePolicy.add("ccdWcsKey", "ccdWcs")
        stagePolicy.add("ampBBoxKey", "ampBBox")

        sourceSet = afwDet.SourceSet()
        for i in xrange(5):
            sourceSet.append(afwDet.Source())
        
        point = afwGeom.makePointD(0.0, 0.0)
        wcs = afwImage.createWcs(point, point, 1, 0, 0, 1);
        ampBBox = afwImage.BBox(afwImage.PointI(0, 0), 1, 1)

        tester = SimpleStageTester(measPipe.SourceToDiaSourceStage(stagePolicy))

        clipboard = pexClipboard.Clipboard()
        clipboard.put(dataPolicy1.get("inputKey"), sourceSet)
        clipboard.put(dataPolicy2.get("inputKey"), sourceSet) 
        clipboard.put(stagePolicy.get("ccdWcsKey"), wcs)
        clipboard.put(stagePolicy.get("ampBBoxKey"), ampBBox)
        
        outWorker = tester.runWorker(clipboard)
    
        for policy in stagePolicy.getPolicyArray("data"):                
            assert(outWorker.contains(policy.get("inputKey")))
            assert(outWorker.contains(policy.get("outputKey")))
            assert(outWorker.contains("persistable_" + policy.get("outputKey")))

            diaSourceSet = outWorker.get(policy.get("outputKey"))
            assert(diaSourceSet.size() == sourceSet.size())

def suite():
    """Returns a suite containing all the test cases in this module."""

    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(SourceToDiaSourceStageTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)

