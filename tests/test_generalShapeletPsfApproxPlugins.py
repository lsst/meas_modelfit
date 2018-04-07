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
import lsst.afw.geom.ellipses
import lsst.afw.table
import lsst.afw.detection
import lsst.log
import lsst.log.utils
import lsst.meas.modelfit
import lsst.meas.base
import lsst.meas.algorithms

#   Set trace to 0-5 to view debug messages.  Level 5 enables all traces.
lsst.log.utils.traceSetAt("meas.modelfit.optimizer.Optimizer", -1)
lsst.log.utils.traceSetAt("meas.modelfit.optimizer.solveTrustRegion", -1)


class GeneralShapeletPsfApproxPluginsTestCase(lsst.utils.tests.TestCase):

    def makeBlankConfig(self):
        config = lsst.meas.base.SingleFrameMeasurementTask.ConfigClass()
        config.slots.centroid = None
        config.slots.shape = None
        config.slots.psfFlux = None
        config.slots.apFlux = None
        config.slots.instFlux = None
        config.slots.modelFlux = None
        config.slots.calibFlux = None
        config.doReplaceWithNoise = False
        return config

    def setUp(self):
        numpy.random.seed(500)
        self.psfSigma = 2.0
        self.exposure = lsst.afw.image.ExposureF(41, 41)
        self.psf = lsst.afw.detection.GaussianPsf(19, 19, self.psfSigma)
        self.schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        self.centroidKey = lsst.afw.table.Point2DKey.addFields(self.schema, "centroid", "centroid", "pixel")
        self.schema.getAliasMap().set("slot_Centroid", "centroid")
        self.psfDir = os.path.join(os.environ["MEAS_MODELFIT_DIR"], "tests", "data", "psfs")

    def tearDown(self):
        del self.exposure
        del self.psf
        del self.schema
        del self.centroidKey
        del self.psfDir

    def checkResult(self, msf):
        # Because we're fitting multiple shapelets to a single Gaussian (a single 0th-order shapelet)
        # we should be able to fit with zero residuals, aside from (single-precision) round-off error.
        dataImage = self.exposure.getPsf().computeImage()
        modelImage = dataImage.Factory(dataImage.getBBox())
        modelImage.getArray()[:, :] *= -1
        msf.evaluate().addToImage(modelImage)
        self.assertFloatsAlmostEqual(dataImage.getArray(), modelImage.getArray(), atol=1E-6,
                                     plotOnFailure=False)

    def testSingleFrame(self):
        self.exposure.setPsf(self.psf)
        config = self.makeBlankConfig()
        config.plugins.names = ["modelfit_GeneralShapeletPsfApprox"]
        config.plugins["modelfit_GeneralShapeletPsfApprox"].sequence = ["SingleGaussian"]
        task = lsst.meas.base.SingleFrameMeasurementTask(config=config, schema=self.schema)
        measCat = lsst.afw.table.SourceCatalog(self.schema)
        measRecord = measCat.addNew()
        measRecord.set(self.centroidKey, lsst.afw.geom.Point2D(20.0, 20.0))
        task.run(measCat, self.exposure)
        keySingleGaussian = lsst.shapelet.MultiShapeletFunctionKey(
            self.schema["modelfit"]["GeneralShapeletPsfApprox"]["SingleGaussian"]
        )
        msfSingleGaussian = measRecord.get(keySingleGaussian)
        self.assertEqual(len(msfSingleGaussian.getComponents()), 1)
        self.checkResult(msfSingleGaussian)

    def testForced(self):
        self.exposure.setPsf(self.psf)
        config = lsst.meas.base.ForcedMeasurementTask.ConfigClass()
        config.slots.centroid = "base_TransformedCentroid"
        config.slots.shape = None
        config.slots.psfFlux = None
        config.slots.apFlux = None
        config.slots.instFlux = None
        config.slots.modelFlux = None
        config.doReplaceWithNoise = False
        config.slots.centroid = "base_TransformedCentroid"
        config.plugins.names = ["base_TransformedCentroid", "modelfit_GeneralShapeletPsfApprox"]
        config.plugins["modelfit_GeneralShapeletPsfApprox"].sequence = ["SingleGaussian"]
        config.copyColumns = {"id": "objectId", "parent": "parentObjectId"}
        refCat = lsst.afw.table.SourceCatalog(self.schema)
        refRecord = refCat.addNew()
        refRecord.set(self.centroidKey, lsst.afw.geom.Point2D(20.0, 20.0))
        refWcs = self.exposure.getWcs()  # same as measurement Wcs
        task = lsst.meas.base.ForcedMeasurementTask(config=config, refSchema=self.schema)
        measCat = task.generateMeasCat(self.exposure, refCat, refWcs)
        task.run(measCat, self.exposure, refCat, refWcs)
        measRecord = measCat[0]
        measSchema = measCat.schema
        keySingleGaussian = lsst.shapelet.MultiShapeletFunctionKey(
            measSchema["modelfit"]["GeneralShapeletPsfApprox"]["SingleGaussian"]
        )
        msfSingleGaussian = measRecord.get(keySingleGaussian)
        self.assertEqual(len(msfSingleGaussian.getComponents()), 1)
        self.checkResult(msfSingleGaussian)

    def testNanFlag(self):
        config = self.makeBlankConfig()
        config.plugins.names = ["modelfit_GeneralShapeletPsfApprox"]
        config.plugins["modelfit_GeneralShapeletPsfApprox"].sequence = ["Full"]
        task = lsst.meas.base.SingleFrameMeasurementTask(config=config, schema=self.schema)
        measCat = lsst.afw.table.SourceCatalog(self.schema)
        measRecord = measCat.addNew()
        psfImage = lsst.afw.image.ImageD(os.path.join(self.psfDir, "galsimPsf_0.9.fits"))
        psfImage.getArray()[0, 0] = numpy.nan
        psfImage.setXY0(lsst.afw.geom.Point2I(0, 0))
        kernel = lsst.afw.math.FixedKernel(psfImage)
        psf = lsst.meas.algorithms.KernelPsf(kernel)
        self.exposure.setPsf(psf)
        center = lsst.afw.geom.Point2D(psfImage.getArray().shape[0]/2, psfImage.getArray().shape[1]/2)
        measRecord.set(self.centroidKey, center)
        task.run(measCat, self.exposure)
        self.assertTrue(measRecord.get("modelfit_GeneralShapeletPsfApprox_Full_flag"))
        self.assertTrue(measRecord.get("modelfit_GeneralShapeletPsfApprox_Full_flag_contains_nan"))
        self.assertFalse(measRecord.get("modelfit_GeneralShapeletPsfApprox_Full_flag_max_inner_iterations"))
        self.assertFalse(measRecord.get("modelfit_GeneralShapeletPsfApprox_Full_flag_max_outer_iterations"))
        self.assertFalse(measRecord.get("modelfit_GeneralShapeletPsfApprox_Full_flag_exception"))

    def testInnerIterationsFlag(self):
        config = self.makeBlankConfig()
        config.plugins.names = ["modelfit_GeneralShapeletPsfApprox"]
        config.plugins["modelfit_GeneralShapeletPsfApprox"].sequence = ["Full"]
        config.plugins["modelfit_GeneralShapeletPsfApprox"].models["Full"].optimizer.maxInnerIterations = 1
        task = lsst.meas.base.SingleFrameMeasurementTask(config=config, schema=self.schema)
        measCat = lsst.afw.table.SourceCatalog(self.schema)
        measRecord = measCat.addNew()
        psfImage = lsst.afw.image.ImageD(os.path.join(self.psfDir, "galsimPsf_0.9.fits"))
        psfImage.setXY0(lsst.afw.geom.Point2I(0, 0))
        kernel = lsst.afw.math.FixedKernel(psfImage)
        psf = lsst.meas.algorithms.KernelPsf(kernel)
        self.exposure.setPsf(psf)
        center = lsst.afw.geom.Point2D(psfImage.getArray().shape[0]/2, psfImage.getArray().shape[1]/2)
        measRecord.set(self.centroidKey, center)
        task.run(measCat, self.exposure)
        self.assertTrue(measRecord.get("modelfit_GeneralShapeletPsfApprox_Full_flag"))
        self.assertFalse(measRecord.get("modelfit_GeneralShapeletPsfApprox_Full_flag_contains_nan"))
        self.assertTrue(measRecord.get("modelfit_GeneralShapeletPsfApprox_Full_flag_max_inner_iterations"))
        self.assertFalse(measRecord.get("modelfit_GeneralShapeletPsfApprox_Full_flag_max_outer_iterations"))
        self.assertFalse(measRecord.get("modelfit_GeneralShapeletPsfApprox_Full_flag_exception"))

    def testOuterIterationsFlag(self):
        config = self.makeBlankConfig()
        config.plugins.names = ["modelfit_GeneralShapeletPsfApprox"]
        config.plugins["modelfit_GeneralShapeletPsfApprox"].sequence = ["Full"]
        config.plugins["modelfit_GeneralShapeletPsfApprox"].models["Full"].optimizer.maxOuterIterations = 1
        task = lsst.meas.base.SingleFrameMeasurementTask(config=config, schema=self.schema)
        measCat = lsst.afw.table.SourceCatalog(self.schema)
        measRecord = measCat.addNew()
        psfImage = lsst.afw.image.ImageD(os.path.join(self.psfDir, "galsimPsf_0.9.fits"))
        psfImage.setXY0(lsst.afw.geom.Point2I(0, 0))
        kernel = lsst.afw.math.FixedKernel(psfImage)
        psf = lsst.meas.algorithms.KernelPsf(kernel)
        self.exposure.setPsf(psf)
        center = lsst.afw.geom.Point2D(psfImage.getArray().shape[0]/2, psfImage.getArray().shape[1]/2)
        measRecord.set(self.centroidKey, center)
        task.run(measCat, self.exposure)
        self.assertTrue(measRecord.get("modelfit_GeneralShapeletPsfApprox_Full_flag"))
        self.assertFalse(measRecord.get("modelfit_GeneralShapeletPsfApprox_Full_flag_contains_nan"))
        self.assertFalse(measRecord.get("modelfit_GeneralShapeletPsfApprox_Full_flag_max_inner_iterations"))
        self.assertTrue(measRecord.get("modelfit_GeneralShapeletPsfApprox_Full_flag_max_outer_iterations"))
        self.assertFalse(measRecord.get("modelfit_GeneralShapeletPsfApprox_Full_flag_exception"))


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
