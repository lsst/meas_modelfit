#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2008-2014 LSST Corporation.
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

import unittest
import os
import numpy
import glob

import lsst.utils.tests
import lsst.shapelet
import lsst.afw.geom.ellipses
import lsst.meas.multifit
import lsst.afw.math
import lsst.afw.image
import lsst.meas.base

import lsst.pex.logging
lsst.pex.logging.Debug("meas.multifit.optimizer", 10)

NUMBER_OF_TEST_KERNELS = 2;
DATA_DIR = os.path.join(os.environ["MEAS_MULTIFIT_DIR"], "tests", "data")

def computeMoments(image):
    """Helper function to compute moments of a postage stamp about its origin."""
    maskedImage = lsst.afw.image.makeMaskedImage(image)
    result = lsst.meas.base.SdssShapeAlgorithm.Result()
    lsst.meas.base.SdssShapeAlgorithm.apply(
        maskedImage,
        lsst.afw.detection.Footprint(image.getBBox(lsst.afw.image.PARENT)),
        lsst.afw.geom.Point2D(0.0, 0.0),
        result
        )
    return result.getShape()

class PsfFitterApplyTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        """Build a dict of configs for testing.
        """
        self.configs = {}
        self.configs['ellipse'] = lsst.meas.multifit.PsfFitterConfig()
        self.configs['ellipse'].primary.ellipticityPriorSigma = 0.3
        self.configs['ellipse'].primary.ellipticityPriorSigma = 0.3
        self.configs['ellipse'].primary.radiusPriorSigma = 0.5
        self.configs['ellipse'].wings.order = 0
        self.configs['ellipse'].wings.ellipticityPriorSigma = 0.0
        self.configs['ellipse'].wings.ellipticityPriorSigma = 0.0
        self.configs['ellipse'].wings.radiusPriorSigma = 0.0

    def tearDown(self):
        del self.configs

    def testApply(self):
        for filename in glob.glob(os.path.join(DATA_DIR, "psfs", "*.fits")):
            kernelImageD = lsst.afw.image.ImageD(filename)
            kernelImageF = kernelImageD.convertF()
            shape = computeMoments(kernelImageF)
            for configKey in self.configs:
                fitter = lsst.meas.multifit.PsfFitter(self.configs[configKey].makeControl())
                multiShapeletFit = fitter.apply(kernelImageF,
                                                0.1,  #not sure what to put here
                                                shape)
                print "Moments:", shape
                for element in multiShapeletFit.getElements():
                    print element.getEllipse().getCore()
                modelImageD = lsst.afw.image.ImageD(kernelImageD.getBBox(lsst.afw.image.PARENT))
                multiShapeletFit.evaluate().addToImage(modelImageD)
                modelImageF = modelImageD.convertF()
                self.assertClose(kernelImageD.getArray(), modelImageD.getArray(), plotOnFailure=True)

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(PsfFitterApplyTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
