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
import unittest

import eups

import lsst.utils.tests as utilsTests
import lsst.daf.base as dafBase
import lsst.pex.policy as pexPolicy
import lsst.afw.detection as afwDet
import lsst.afw.image as afwImage
import lsst.meas.pipeline as measPipe

from lsst.pex.harness.Clipboard import Clipboard
from lsst.pex.harness.simpleStageTester import SimpleStageTester

class ComputeSourceSkyCoordsStageTestCase(unittest.TestCase):
    """Tests the ComputeSourceSkyCoordsStage pipeline stage.
    """

    def setUp(self):
        self.sources = afwDet.SourceSet()
        for i in xrange(1, 10):
            s = afwDet.Source()
            s.setXFlux(i*10.0)
            s.setXFluxErr(i * 0.01)
            s.setYFlux(i*10.0)
            s.setYFluxErr(i * 0.01)
            s.setXAstrom(i*10.0)
            s.setXAstromErr(i * 0.01)
            s.setYAstrom(i*10.0)
            s.setYAstromErr(i * 0.01)
            s.setXPeak(i*10.0)
            s.setYPeak(i*10.0)
            self.sources.append(s)
        imagePath = os.path.join(eups.productDir("meas_pipeline"), "tests",
                                 "v695833-e0-c000-a00.sci_img.fits")
        self.wcs = afwImage.makeWcs(afwImage.readMetadata(imagePath))

    def tearDown(self):
        del self.sources
        del self.wcs

    def testStage(self):
        # setup stage environment
        policyFile = pexPolicy.DefaultPolicyFile(
            "meas_pipeline", "ComputeSourceSkyCoordsStageDictionary.paf", "policy")
        policy = pexPolicy.Policy.createPolicy(
            policyFile, policyFile.getRepositoryPath())
        clipboard = Clipboard()
        clipboard.put(policy.get("inputKeys.wcs"), self.wcs)
        clipboard.put(policy.get("inputKeys.sources"), self.sources)
        stage = measPipe.ComputeSourceSkyCoordsStage(policy)

        # run stage
        tester = SimpleStageTester(stage)
        output = tester.runWorker(clipboard)

        # verify output
        self.assertTrue(output.contains(policy.get("inputKeys.wcs")))
        self.assertTrue(output.contains(policy.get("inputKeys.sources")))
        sources = output.get(policy.get("inputKeys.sources"))
        self.assertTrue(isinstance(sources, afwDet.SourceSet))
        for s in sources:
            self.assertNotEqual(s.getRaFlux(), 0.0)
            self.assertNotEqual(s.getRaFluxErr(), 0.0)
            self.assertNotEqual(s.getDecFlux(), 0.0)
            self.assertNotEqual(s.getDecFluxErr(), 0.0)
            self.assertNotEqual(s.getRaAstrom(), 0.0)
            self.assertNotEqual(s.getRaAstromErr(), 0.0)
            self.assertNotEqual(s.getDecAstrom(), 0.0)
            self.assertNotEqual(s.getDecAstromErr(), 0.0)
            self.assertNotEqual(s.getRaPeak(), 0.0)
            self.assertNotEqual(s.getDecPeak(), 0.0)
            self.assertEqual(s.getRa(), s.getRaAstrom())
            self.assertEqual(s.getRaErrForDetection(), s.getRaAstromErr())
            self.assertEqual(s.getDec(), s.getDecAstrom())
            self.assertEqual(s.getDecErrForDetection(), s.getDecAstromErr())


def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()
    suites = [unittest.makeSuite(ComputeSourceSkyCoordsStageTestCase),
              unittest.makeSuite(utilsTests.MemoryTestCase),
             ]
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == '__main__':
    run(True)
 
