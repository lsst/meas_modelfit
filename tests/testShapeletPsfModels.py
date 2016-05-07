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
import numpy
import os
import time

import lsst.utils.tests
import lsst.shapelet
import lsst.afw.geom.ellipses
import lsst.afw.table
import lsst.afw.detection
import lsst.meas.modelfit
import lsst.meas.base
import lsst.meas.algorithms

numpy.random.seed(500)

lsst.pex.logging.Debug("meas.modelfit.optimizer.Optimizer", 0)
lsst.pex.logging.Debug("meas.modelfit.optimizer.solveTrustRegion", 0)

class ShapeletPsfApproxPluginsTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.exposure = lsst.afw.image.ExposureF(41, 41)
        self.schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        self.centroidKey = lsst.afw.table.Point2DKey.addFields(self.schema, "centroid", "centroid", "pixel")
        self.schema.getAliasMap().set("slot_Centroid", "centroid")
        self.psfDir = os.path.join(os.environ["MEAS_MODELFIT_DIR"], "tests", "data", "psfs")

    def tearDown(self):
        del self.exposure
        del self.schema
        del self.centroidKey
        del self.psfDir

    def makePsf(self, psfname, max=None):
        data = lsst.afw.image.ImageF(os.path.join(self.psfDir, psfname)).getArray().astype(numpy.float64)
        if max is not None:
            trim0 = (data.shape[0] - max)/2
            trim1 = (data.shape[1] - max)/2
            if trim0 > 0 and trim1 > 0:
                data = data[trim0:trim0+max, trim1:trim1+max]
        kernel = lsst.afw.math.FixedKernel(lsst.afw.image.ImageD(data))
        return lsst.meas.algorithms.KernelPsf(kernel)

    def runTask(self, psftype, sequence):
        config = lsst.meas.base.SingleFrameMeasurementTask.ConfigClass()
        config.slots.centroid = None
        config.slots.shape = None
        config.slots.psfFlux = None
        config.slots.apFlux = None
        config.slots.instFlux = None
        config.slots.modelFlux = None
        config.slots.calibFlux = None
        config.doReplaceWithNoise = False
        config.plugins.names = ["modelfit_ShapeletPsfApprox"]
        config.plugins["modelfit_ShapeletPsfApprox"].sequence = sequence
        task = lsst.meas.base.SingleFrameMeasurementTask(config=config, schema=self.schema)
        measCat = lsst.afw.table.SourceCatalog(self.schema)
        measRecord = measCat.addNew()
        measRecord.set(self.centroidKey, lsst.afw.geom.Point2D(20.0, 20.0))
        measRecord.set(self.centroidKey, lsst.afw.geom.Point2D(20.0, 20.0))
        task.run(measCat, self.exposure)
        return measRecord

    def testSingleGaussian(self):
        sigma1 = 3.0
        self.exposure.setPsf(lsst.afw.detection.GaussianPsf(19, 19, sigma1))
        measRecord = self.runTask("Single Gaussian Psf", ["SingleGaussian"])
        keySingleGaussian = lsst.shapelet.MultiShapeletFunctionKey(
            self.schema["modelfit"]["ShapeletPsfApprox"]["SingleGaussian"]
            )
        msfSingleGaussian = measRecord.get(keySingleGaussian)
        self.assertEqual(len(msfSingleGaussian.getComponents()), 1)
        comps = msfSingleGaussian.getComponents()
        r0 = comps[0].getEllipse().getCore().getDeterminantRadius()
        self.assertClose(r0, sigma1, .05)

    def testDoubleGaussian(self):
        sigma1 = 2.0
        sigma2 = 4.0
        self.exposure.setPsf(lsst.meas.algorithms.DoubleGaussianPsf(19, 19, sigma1, sigma2, .25))
        measRecord = self.runTask("Double Gaussian Psf", ["DoubleGaussian"])
        keyDoubleGaussian = lsst.shapelet.MultiShapeletFunctionKey(
            self.schema["modelfit"]["ShapeletPsfApprox"]["DoubleGaussian"]
            )
        msf = measRecord.get(keyDoubleGaussian)
        comps = msf.getComponents()
        self.assertEqual(len(comps), 2)
        #  amplitudes are equal by construction
        A0 = measRecord.get("modelfit_ShapeletPsfApprox_DoubleGaussian_0_0")
        A1 = measRecord.get("modelfit_ShapeletPsfApprox_DoubleGaussian_1_0")
        self.assertClose(A0, A1, .05)
        r0 = comps[0].getEllipse().getCore().getDeterminantRadius()
        r1 = comps[1].getEllipse().getCore().getDeterminantRadius()
        self.assertClose(r0, sigma1, .05)
        self.assertClose(r1, sigma2, .05)

    def testDoubleShapelet(self):
        self.exposure.setPsf(self.makePsf("galsimPsf_0.5.fits", max=33))
        measRecord = self.runTask("Galsim Psf", ["DoubleShapelet"])
        keyDoubleShapelet = lsst.shapelet.MultiShapeletFunctionKey(
            self.schema["modelfit"]["ShapeletPsfApprox"]["DoubleShapelet"]
            )
        msf = measRecord.get(keyDoubleShapelet)
        comps = msf.getComponents()
        self.assertEqual(len(comps), 2)
        A0 = measRecord.get("modelfit_ShapeletPsfApprox_DoubleShapelet_0_0")
        A1 = measRecord.get("modelfit_ShapeletPsfApprox_DoubleShapelet_1_0")
        self.assertGreater(A1, .04)
        self.assertGreater(A0, .04)

    def testFull(self):
        self.exposure.setPsf(self.makePsf("galsimPsf_0.9.fits", max=33))
        measRecord = self.runTask("Galsim Psf", ["Full"])
        keyFull = lsst.shapelet.MultiShapeletFunctionKey(
            self.schema["modelfit"]["ShapeletPsfApprox"]["Full"]
            )
        msf = measRecord.get(keyFull)
        comps = msf.getComponents()
        self.assertEqual(len(comps), 4)
        A1 = measRecord.get("modelfit_ShapeletPsfApprox_Full_1_0")
        A2 = measRecord.get("modelfit_ShapeletPsfApprox_Full_2_0")
        # test the primary and wings to be sure we are getting something
        self.assertGreater(A2, .04)
        self.assertGreater(A1, .04)

    def testSequence(self):
        sigma1 = 2.0
        sigma2 = 4.0
        self.exposure.setPsf(lsst.meas.algorithms.DoubleGaussianPsf(19, 19, sigma1, sigma2, .25))
        measRecord = self.runTask("Single Gaussian Psf", ["SingleGaussian", "DoubleGaussian",
                                  "DoubleShapelet"])
        keySingleGaussian = lsst.shapelet.MultiShapeletFunctionKey(
            self.schema["modelfit"]["ShapeletPsfApprox"]["SingleGaussian"]
            )
        msfSingleGaussian = measRecord.get(keySingleGaussian)
        self.assertEqual(len(msfSingleGaussian.getComponents()), 1)
        comps = msfSingleGaussian.getComponents()
        r0 = comps[0].getEllipse().getCore().getDeterminantRadius()
        # don't expect it to be all that close, but the DoubleGaussian should be
        self.assertClose(r0, sigma1, .3)

        keyDoubleGaussian = lsst.shapelet.MultiShapeletFunctionKey(
            self.schema["modelfit"]["ShapeletPsfApprox"]["DoubleGaussian"]
            )
        msfDoubleGaussian = measRecord.get(keyDoubleGaussian)
        comps = msfDoubleGaussian.getComponents()
        r0 = comps[0].getEllipse().getCore().getDeterminantRadius()
        r1 = comps[1].getEllipse().getCore().getDeterminantRadius()
        self.assertClose(r0, sigma1, .05)
        self.assertClose(r1, sigma2, .05)

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(ShapeletPsfApproxPluginsTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
