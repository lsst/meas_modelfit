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
import lsst.utils.tests
import lsst.shapelet
import lsst.afw.geom.ellipses
import lsst.meas.multifit
import lsst.afw.math
import lsst.afw.image

NUMBER_OF_TEST_KERNELS = 2;
DATA_DIR = os.path.join(os.environ["MEAS_MULTIFIT_DIR"], "tests", "data")

class PsfFitterApplyTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        """Build a dict of configs for testing.
        This was directly cut and pasted from testPsfFitter. 
        I will either extend these here, or add the testApply to testPsfFitter.py.
        """
        self.configs = {}
        self.configs['fixed'] = lsst.meas.multifit.PsfFitterConfig()
        self.configs['ellipse'] = lsst.meas.multifit.PsfFitterConfig()
        self.configs['ellipse'].primary.ellipticityPriorSigma = 0.3
        self.configs['ellipse'].primary.ellipticityPriorSigma = 0.3
        self.configs['ellipse'].primary.radiusPriorSigma = 0.5
        self.configs['ellipse'].wings.ellipticityPriorSigma = 0.3
        self.configs['ellipse'].wings.ellipticityPriorSigma = 0.3
        self.configs['ellipse'].wings.radiusPriorSigma = 0.5
        self.configs['full'] = lsst.meas.multifit.PsfFitterConfig()
        self.configs['full'].inner.order = 0
        self.configs['full'].primary.order = 4
        self.configs['full'].wings.order = 4
        self.configs['full'].outer.order = 0
        self.configs['full'].inner.ellipticityPriorSigma = 0.3
        self.configs['full'].inner.ellipticityPriorSigma = 0.3
        self.configs['full'].inner.radiusPriorSigma = 0.5
        self.configs['full'].inner.positionPriorSigma = 0.1
        self.configs['full'].primary.ellipticityPriorSigma = 0.3
        self.configs['full'].primary.ellipticityPriorSigma = 0.3
        self.configs['full'].primary.radiusPriorSigma = 0.5
        self.configs['full'].primary.positionPriorSigma = 0.1
        self.configs['full'].wings.ellipticityPriorSigma = 0.3
        self.configs['full'].wings.ellipticityPriorSigma = 0.3
        self.configs['full'].wings.radiusPriorSigma = 0.5
        self.configs['full'].wings.positionPriorSigma = 0.1
        self.configs['full'].outer.ellipticityPriorSigma = 0.3
        self.configs['full'].outer.ellipticityPriorSigma = 0.3
        self.configs['full'].outer.radiusPriorSigma = 0.5
        self.configs['full'].outer.positionPriorSigma = 0.1
        self.sctrl = lsst.afw.math.StatisticsControl()

    def tearDown(self):
        del self.configs

    def testApply(self):
        for idxKernel in range(NUMBER_OF_TEST_KERNELS):
            kernel = lsst.meas.algorithms.KernelPsfPersistableFacade_readFits(
                         os.path.join(DATA_DIR, 'psfKernel%i.fits'%(idxKernel)))
            kernelImage = kernel.computeKernelImage().convertF()
            for configKey in self.configs:
                fitter = lsst.meas.multifit.PsfFitter(self.configs[configKey].makeControl())
                multiShapeletFit = fitter.apply(kernelImage,
                                                0.05,  #not sure what to put here
                                                kernel.computeShape())
                residualImage = -lsst.afw.image.ImageF(kernelImage)
                multiShapeletFit.evaluate().addToImage(residualImage)
                self.assertClose(0.0, lsst.afw.math.makeStatistics(residualImage,
                                                                   lsst.afw.math.MEANSQUARE,
                                                                   self.sctr).getValue())

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
