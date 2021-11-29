#
# LSST Data Management System
#
# Copyright 2008-2016  AURA/LSST.
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
# see <https://www.lsstcorp.org/LegalNotices/>.
#
import unittest
import numpy
import os

import lsst.utils.tests
import lsst.shapelet
import lsst.geom
import lsst.afw.geom.ellipses
import lsst.afw.table
import lsst.afw.detection
import lsst.log
import lsst.utils.logging
import lsst.meas.modelfit
import lsst.meas.base
import lsst.meas.algorithms

#   Set trace to 0-5 to view debug messages.  Level 5 enables all traces.
lsst.utils.logging.trace_set_at("meas.modelfit.optimizer.Optimizer", -1)
lsst.utils.logging.trace_set_at("meas.modelfit.optimizer.solveTrustRegion", -1)


class GeneralShapeletPsfApproxPluginsTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        numpy.random.seed(500)
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
            trim0 = (data.shape[0] - max)//2
            trim1 = (data.shape[1] - max)//2
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
        config.slots.gaussianFlux = None
        config.slots.modelFlux = None
        config.slots.calibFlux = None
        config.doReplaceWithNoise = False
        config.plugins.names = ["modelfit_GeneralShapeletPsfApprox"]
        config.plugins["modelfit_GeneralShapeletPsfApprox"].sequence = sequence
        task = lsst.meas.base.SingleFrameMeasurementTask(config=config, schema=self.schema)
        measCat = lsst.afw.table.SourceCatalog(self.schema)
        measRecord = measCat.addNew()
        measRecord.set(self.centroidKey, lsst.geom.Point2D(20.0, 20.0))
        measRecord.set(self.centroidKey, lsst.geom.Point2D(20.0, 20.0))
        task.run(measCat, self.exposure)
        return measRecord

    def testSingleGaussian(self):
        sigma1 = 3.0
        self.exposure.setPsf(lsst.afw.detection.GaussianPsf(19, 19, sigma1))
        measRecord = self.runTask("Single Gaussian Psf", ["SingleGaussian"])
        keySingleGaussian = lsst.shapelet.MultiShapeletFunctionKey(
            self.schema["modelfit"]["GeneralShapeletPsfApprox"]["SingleGaussian"]
        )
        msfSingleGaussian = measRecord.get(keySingleGaussian)
        self.assertEqual(len(msfSingleGaussian.getComponents()), 1)
        comps = msfSingleGaussian.getComponents()
        r0 = comps[0].getEllipse().getCore().getDeterminantRadius()
        self.assertFloatsAlmostEqual(r0, sigma1, .05)

    def testDoubleGaussian(self):
        sigma1 = 2.0
        sigma2 = 4.0
        self.exposure.setPsf(lsst.meas.algorithms.DoubleGaussianPsf(19, 19, sigma1, sigma2, .25))
        measRecord = self.runTask("Double Gaussian Psf", ["DoubleGaussian"])
        keyDoubleGaussian = lsst.shapelet.MultiShapeletFunctionKey(
            self.schema["modelfit"]["GeneralShapeletPsfApprox"]["DoubleGaussian"]
        )
        msf = measRecord.get(keyDoubleGaussian)
        comps = msf.getComponents()
        self.assertEqual(len(comps), 2)
        #  amplitudes are equal by construction
        A0 = measRecord.get("modelfit_GeneralShapeletPsfApprox_DoubleGaussian_0_0")
        A1 = measRecord.get("modelfit_GeneralShapeletPsfApprox_DoubleGaussian_1_0")
        self.assertFloatsAlmostEqual(A0, A1, .05)
        r0 = comps[0].getEllipse().getCore().getDeterminantRadius()
        r1 = comps[1].getEllipse().getCore().getDeterminantRadius()
        self.assertFloatsAlmostEqual(r0, sigma1, .05)
        self.assertFloatsAlmostEqual(r1, sigma2, .05)

    def testDoubleShapelet(self):
        self.exposure.setPsf(self.makePsf("galsimPsf_0.5.fits", max=33))
        measRecord = self.runTask("Galsim Psf", ["DoubleShapelet"])
        keyDoubleShapelet = lsst.shapelet.MultiShapeletFunctionKey(
            self.schema["modelfit"]["GeneralShapeletPsfApprox"]["DoubleShapelet"]
        )
        msf = measRecord.get(keyDoubleShapelet)
        comps = msf.getComponents()
        self.assertEqual(len(comps), 2)
        A0 = measRecord.get("modelfit_GeneralShapeletPsfApprox_DoubleShapelet_0_0")
        A1 = measRecord.get("modelfit_GeneralShapeletPsfApprox_DoubleShapelet_1_0")
        self.assertGreater(A1, .04)
        self.assertGreater(A0, .04)

    def testFull(self):
        self.exposure.setPsf(self.makePsf("galsimPsf_0.9.fits", max=33))
        measRecord = self.runTask("Galsim Psf", ["Full"])
        keyFull = lsst.shapelet.MultiShapeletFunctionKey(
            self.schema["modelfit"]["GeneralShapeletPsfApprox"]["Full"]
        )
        msf = measRecord.get(keyFull)
        comps = msf.getComponents()
        self.assertEqual(len(comps), 4)
        A1 = measRecord.get("modelfit_GeneralShapeletPsfApprox_Full_1_0")
        A2 = measRecord.get("modelfit_GeneralShapeletPsfApprox_Full_2_0")
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
            self.schema["modelfit"]["GeneralShapeletPsfApprox"]["SingleGaussian"]
        )
        msfSingleGaussian = measRecord.get(keySingleGaussian)
        self.assertEqual(len(msfSingleGaussian.getComponents()), 1)
        comps = msfSingleGaussian.getComponents()
        r0 = comps[0].getEllipse().getCore().getDeterminantRadius()
        # don't expect it to be all that close, but the DoubleGaussian should be
        self.assertFloatsAlmostEqual(r0, sigma1, .3)

        keyDoubleGaussian = lsst.shapelet.MultiShapeletFunctionKey(
            self.schema["modelfit"]["GeneralShapeletPsfApprox"]["DoubleGaussian"]
        )
        msfDoubleGaussian = measRecord.get(keyDoubleGaussian)
        comps = msfDoubleGaussian.getComponents()
        r0 = comps[0].getEllipse().getCore().getDeterminantRadius()
        r1 = comps[1].getEllipse().getCore().getDeterminantRadius()
        self.assertFloatsAlmostEqual(r0, sigma1, .05)
        self.assertFloatsAlmostEqual(r1, sigma2, .05)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
