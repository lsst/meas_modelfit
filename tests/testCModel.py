#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2008-2013 LSST Corporation.
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
import lsst.afw.image
import lsst.meas.modelfit
import lsst.meas.base

numpy.random.seed(500)

lsst.pex.logging.Debug("meas.modelfit.optimizer.Optimizer", 0)
lsst.pex.logging.Debug("meas.modelfit.optimizer.solveTrustRegion", 0)

def makeMultiShapeletCircularGaussian(sigma):
    s = lsst.shapelet.ShapeletFunction(0, lsst.shapelet.HERMITE, sigma)
    s.getCoefficients()[0] = 1.0 / lsst.shapelet.ShapeletFunction.FLUX_FACTOR
    m = lsst.shapelet.MultiShapeletFunction()
    m.getComponents().push_back(s)
    return m

def computePsfFlux(centroid, exposure):
    schema = lsst.afw.table.SourceTable.makeMinimalSchema()
    pointKey = lsst.afw.table.Point2DKey.addFields(schema, "centroid", "known input centroid", "pixels")
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
        crval = lsst.afw.coord.IcrsCoord(45.0*lsst.afw.geom.degrees, 45.0*lsst.afw.geom.degrees)
        crpix = lsst.afw.geom.Point2D(0.0, 0.0)
        cdelt = (0.2*lsst.afw.geom.arcseconds).asDegrees()
        dataWcs = lsst.afw.image.makeWcs(crval, crpix, cdelt, 0.0, 0.0, cdelt)
        dataCalib = lsst.afw.image.Calib()
        dataCalib.setFluxMag0(1e12)
        self.xyPosition = lsst.afw.geom.Point2D(1.1, -0.8)
        position = dataWcs.pixelToSky(self.xyPosition)
        bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-100, -100), lsst.afw.geom.Point2I(100, 100))
        self.exposure = lsst.afw.image.ExposureF(bbox)
        self.exposure.setWcs(dataWcs)
        self.exposure.setCalib(dataCalib)
        self.trueFlux = 65.0
        self.psfSigma = 2.0
        psf = lsst.afw.detection.GaussianPsf(25, 25, self.psfSigma)
        self.exposure.setPsf(psf)
        psfImage = psf.computeImage(self.xyPosition)
        psfImage.getArray()[:,:] *= self.trueFlux
        psfBBox = psfImage.getBBox(lsst.afw.image.PARENT)
        self.footprint = lsst.afw.detection.Footprint(psfBBox)
        subImage = lsst.afw.image.ImageF(self.exposure.getMaskedImage().getImage(), psfBBox,
                                         lsst.afw.image.PARENT)
        subImage.getArray()[:,:] = psfImage.getArray()

    def tearDown(self):
        del self.xyPosition
        del self.exposure
        del self.trueFlux
        del self.psfSigma
        del self.footprint

    def testNoNoise(self):
        """Test that CModelAlgorithm.apply() works when applied to a postage-stamp
        containing only a point source with no noise.
        """
        ctrl = lsst.meas.modelfit.CModelControl()
        algorithm = lsst.meas.modelfit.CModelAlgorithm(ctrl)
        result = algorithm.apply(
            self.exposure, makeMultiShapeletCircularGaussian(self.psfSigma),
            self.xyPosition, self.exposure.getPsf().computeShape()
            )
        self.assertFalse(result.initial.getFlag(result.FAILED))
        self.assertClose(result.initial.flux, self.trueFlux, rtol=0.01)
        self.assertClose(result.initial.fluxSigma, 0.0, rtol=0.0)
        self.assertLess(result.initial.getEllipse().getDeterminantRadius(), 0.2)
        self.assertFalse(result.exp.getFlag(result.FAILED))
        self.assertClose(result.exp.flux, self.trueFlux, rtol=0.01)
        self.assertClose(result.exp.fluxSigma, 0.0, rtol=0.0)
        self.assertLess(result.exp.getEllipse().getDeterminantRadius(), 0.2)
        self.assertFalse(result.dev.getFlag(result.FAILED))
        self.assertClose(result.dev.flux, self.trueFlux, rtol=0.01)
        self.assertClose(result.dev.fluxSigma, 0.0, rtol=0.0)
        self.assertLess(result.dev.getEllipse().getDeterminantRadius(), 0.2)
        self.assertFalse(result.getFlag(result.FAILED))
        self.assertClose(result.flux, self.trueFlux, rtol=0.01)
        self.assertClose(result.fluxSigma, 0.0, rtol=0.0, atol=0.0)

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
            self.assertClose(psfFlux, cmodel.flux, rtol=0.1/fluxFactor**0.5)
            self.assertClose(psfFluxSigma, cmodel.fluxSigma, rtol=0.1/fluxFactor**0.5)


def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(CModelTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
