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

import lsst.utils.tests
import lsst.shapelet
import lsst.afw.geom
import lsst.afw.image
import lsst.log
import lsst.log.utils
import lsst.meas.modelfit
import lsst.meas.base

#   Set trace to 0-5 to view debug messages.  Level 5 enables all traces.
lsst.log.utils.traceSetAt("meas.modelfit.optimizer.Optimizer", -1)
lsst.log.utils.traceSetAt("meas.modelfit.optimizer.solveTrustRegion", -1)


def makeMultiShapeletCircularGaussian(sigma):
    s = lsst.shapelet.ShapeletFunction(0, lsst.shapelet.HERMITE, sigma)
    s.getCoefficients()[0] = 1.0 / lsst.shapelet.ShapeletFunction.FLUX_FACTOR
    m = lsst.shapelet.MultiShapeletFunction()
    m.addComponent(s)
    return m


def computePsfFlux(centroid, exposure):
    schema = lsst.afw.table.SourceTable.makeMinimalSchema()
    pointKey = lsst.afw.table.Point2DKey.addFields(schema, "centroid", "known input centroid", "pixel")
    schema.getAliasMap().set("slot_Centroid", "centroid")
    algorithm = lsst.meas.base.PsfFluxAlgorithm(lsst.meas.base.PsfFluxControl(), "base_PsfFlux", schema)
    table = lsst.afw.table.SourceTable.make(schema)
    record = table.makeRecord()
    record.set(pointKey, centroid)
    algorithm.measure(record, exposure)
    return record.get("base_PsfFlux_flux"), record.get("base_PsfFlux_fluxSigma")


class CModelTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        # Setup test data: a single point source, initially with no noise.
        numpy.random.seed(500)
        crval = lsst.afw.geom.SpherePoint(45.0, 45.0, lsst.afw.geom.degrees)
        crpix = lsst.afw.geom.Point2D(0.0, 0.0)
        scale = 0.2 * lsst.afw.geom.arcseconds
        cdMatrix = lsst.afw.geom.makeCdMatrix(scale=scale, flipX=True)
        dataWcs = lsst.afw.geom.makeSkyWcs(crpix=crpix, crval=crval, cdMatrix=cdMatrix)
        dataCalib = lsst.afw.image.Calib()
        dataCalib.setFluxMag0(1e12)
        self.xyPosition = lsst.afw.geom.Point2D(1.1, -0.8)
        bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-100, -100), lsst.afw.geom.Point2I(100, 100))
        self.exposure = lsst.afw.image.ExposureF(bbox)
        self.exposure.setWcs(dataWcs)
        self.exposure.setCalib(dataCalib)
        self.trueFlux = 65.0
        self.psfSigma = 2.0
        psf = lsst.afw.detection.GaussianPsf(25, 25, self.psfSigma)
        self.exposure.setPsf(psf)
        psfImage = psf.computeImage(self.xyPosition)
        psfImage.getArray()[:, :] *= self.trueFlux
        psfBBox = psfImage.getBBox(lsst.afw.image.PARENT)
        subImage = lsst.afw.image.ImageF(self.exposure.getMaskedImage().getImage(), psfBBox,
                                         lsst.afw.image.PARENT)
        subImage.getArray()[:, :] = psfImage.getArray()

    def tearDown(self):
        del self.xyPosition
        del self.exposure
        del self.trueFlux
        del self.psfSigma

    def testNoNoise(self):
        """Test that CModelAlgorithm.apply() works when applied to a postage-stamp
        containing only a point source with no noise.

        We still have to pretend there is noise (i.e. have nonzero values in
        the variance plane) to allow it to compute a likelihood, though.
        """
        ctrl = lsst.meas.modelfit.CModelControl()
        ctrl.initial.usePixelWeights = False
        algorithm = lsst.meas.modelfit.CModelAlgorithm(ctrl)
        var = 1E-16
        self.exposure.getMaskedImage().getVariance().getArray()[:, :] = var
        psfImage = self.exposure.getPsf().computeKernelImage(self.xyPosition).getArray()
        expectedFluxSigma = var**0.5 * (psfImage**2).sum()**(-0.5)
        result = algorithm.apply(
            self.exposure, makeMultiShapeletCircularGaussian(self.psfSigma),
            self.xyPosition, self.exposure.getPsf().computeShape()
        )
        self.assertFalse(result.initial.flags[result.FAILED])
        self.assertFloatsAlmostEqual(result.initial.flux, self.trueFlux, rtol=0.01)
        self.assertFloatsAlmostEqual(result.initial.fluxSigma, expectedFluxSigma, rtol=0.01)
        self.assertLess(result.initial.ellipse.getDeterminantRadius(), 0.2)
        self.assertFalse(result.exp.flags[result.FAILED])
        self.assertFloatsAlmostEqual(result.exp.flux, self.trueFlux, rtol=0.01)
        self.assertFloatsAlmostEqual(result.exp.fluxSigma, expectedFluxSigma, rtol=0.01)
        self.assertLess(result.exp.ellipse.getDeterminantRadius(), 0.2)
        self.assertFalse(result.dev.flags[result.FAILED])
        self.assertFloatsAlmostEqual(result.dev.flux, self.trueFlux, rtol=0.01)
        self.assertFloatsAlmostEqual(result.dev.fluxSigma, expectedFluxSigma, rtol=0.01)
        self.assertLess(result.dev.ellipse.getDeterminantRadius(), 0.2)
        self.assertFalse(result.flags[result.FAILED])
        self.assertFloatsAlmostEqual(result.flux, self.trueFlux, rtol=0.01)

    def testVsPsfFlux(self):
        """Test that CModel produces results comparable to PsfFlux when run
        on point sources.
        """
        noiseSigma = 1.0
        for fluxFactor in (1.0, 10.0, 100.0):
            exposure = self.exposure.Factory(self.exposure, True)
            exposure.getMaskedImage().getImage().getArray()[:] *= fluxFactor
            exposure.getMaskedImage().getVariance().getArray()[:] = noiseSigma**2
            exposure.getMaskedImage().getImage().getArray()[:] += \
                noiseSigma*numpy.random.randn(exposure.getHeight(), exposure.getWidth())
            ctrl = lsst.meas.modelfit.CModelControl()
            algorithm = lsst.meas.modelfit.CModelAlgorithm(ctrl)
            cmodel = algorithm.apply(
                exposure, makeMultiShapeletCircularGaussian(self.psfSigma),
                self.xyPosition, self.exposure.getPsf().computeShape()
            )
            psfFlux, psfFluxSigma = computePsfFlux(self.xyPosition, exposure)
            self.assertFloatsAlmostEqual(psfFlux, cmodel.flux, rtol=0.1/fluxFactor**0.5)
            self.assertFloatsAlmostEqual(psfFluxSigma, cmodel.fluxSigma, rtol=0.1/fluxFactor**0.5)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
