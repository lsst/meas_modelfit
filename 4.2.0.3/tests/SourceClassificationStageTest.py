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

import os, os.path
import pdb
import unittest

import eups
import lsst.utils.tests as utilsTests

from lsst.pex.harness.Clipboard import  Clipboard
import lsst.pex.harness.stage as harnessStage
from lsst.pex.harness.simpleStageTester import SimpleStageTester

import lsst.pex.policy as pexPolicy
import lsst.afw.detection as afwDet
import lsst.meas.pipeline as measPipe
from lsst.meas.utils.sourceClassifiers import *

class SourceClassificationStageTestCase(unittest.TestCase):
    """A test case for the SourceClassificationStage"""

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def testStage(self):
        file = pexPolicy.DefaultPolicyFile("meas_pipeline", 
                "sourceClassification0_policy.paf", "tests")
        policy = pexPolicy.Policy.createPolicy(file)
    
        tester = SimpleStageTester(measPipe.SourceClassificationStage(policy))

        set0 = afwDet.DiaSourceSet()
        set1 = afwDet.DiaSourceSet()
        set0.append(afwDet.DiaSource())
        set0.append(afwDet.DiaSource())
        set0.append(afwDet.DiaSource())
        set1.append(afwDet.DiaSource())
        set1.append(afwDet.DiaSource())
        set1.append(afwDet.DiaSource())

        clipboard = Clipboard()
        clipboard.put("sourceSet0", set0)
        clipboard.put("sourceSet1", set1)

        outWorker = tester.runWorker(clipboard)

        #TODO: insert assert statements here!

def suite():
    """Returns a suite containing all the test cases in this module."""

    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(SourceClassificationStageTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)

