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
import lsst.meas.multifit

numpy.random.seed(500)

lsst.pex.logging.Debug("meas.multifit.optimizer.Optimizer", 0)
lsst.pex.logging.Debug("meas.multifit.optimizer.solveTrustRegion", 0)

def makeMultiShapeletCircularGaussian(sigma):
    s = lsst.shapelet.ShapeletFunction(0, lsst.shapelet.HERMITE, sigma)
    s.getCoefficients()[0] = 1.0 / lsst.shapelet.ShapeletFunction.FLUX_FACTOR
    m = lsst.shapelet.MultiShapeletFunction()
    m.getElements().push_back(s)
    return m


class CModelTestCase(lsst.utils.tests.TestCase):

    def testPointSource(self):
        """Test that CModelAlgorithm.apply() works when applied to a postage-stamp
        containing only a point source
        """
        crval = lsst.afw.coord.IcrsCoord(45.0*lsst.afw.geom.degrees, 45.0*lsst.afw.geom.degrees)
        crpix = lsst.afw.geom.Point2D(0.0, 0.0)
        cdelt = (0.2*lsst.afw.geom.arcseconds).asDegrees()
        dataWcs = lsst.afw.image.makeWcs(crval, crpix, cdelt, 0.0, 0.0, cdelt)
        dataCalib = lsst.afw.image.Calib()
        dataCalib.setFluxMag0(1e12)
        xyPosition = lsst.afw.geom.Point2D(1.1, -0.8)
        position = dataWcs.pixelToSky(xyPosition)
        noiseSigma = 0.01
        bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-100, -100), lsst.afw.geom.Point2I(100, 100))
        exposure = lsst.afw.image.ExposureF(bbox)
        exposure.setWcs(dataWcs)
        exposure.setCalib(dataCalib)
        exposure.getMaskedImage().getVariance().getArray()[:,:] = noiseSigma**2

        # Start by creating a point source and fitting it with an extended model
        trueFlux = 65.0
        exposure.getMaskedImage().getImage().getArray()[:,:] \
            = numpy.random.randn(bbox.getHeight(), bbox.getWidth()) * noiseSigma
        psfSigma = 2.0
        psf = lsst.afw.detection.GaussianPsf(25, 25, psfSigma)
        exposure.setPsf(psf)
        psfImage = psf.computeImage(xyPosition)
        psfImage.getArray()[:,:] *= trueFlux
        psfBBox = psfImage.getBBox(lsst.afw.image.PARENT)
        subImage = lsst.afw.image.ImageF(exposure.getMaskedImage().getImage(), psfBBox, lsst.afw.image.PARENT)
        subImage.getArray()[:,:] = psfImage.getArray()

        ctrl = lsst.meas.multifit.CModelControl()
        algorithm = lsst.meas.multifit.CModelAlgorithm(ctrl)
        result = algorithm.apply(
            exposure, lsst.afw.detection.Footprint(psfBBox),
            makeMultiShapeletCircularGaussian(psfSigma),
            xyPosition, psf.computeShape()
            )
        self.assertFalse(result.initial.getFlag(result.FAILED))
        self.assertClose(result.initial.flux, trueFlux, rtol=0.01)
        self.assertLess(result.initial.getEllipse().getDeterminantRadius(), 0.2)
        self.assertFalse(result.exp.getFlag(result.FAILED))
        self.assertClose(result.exp.flux, trueFlux, rtol=0.01)
        self.assertLess(result.exp.getEllipse().getDeterminantRadius(), 0.2)
        self.assertFalse(result.dev.getFlag(result.FAILED))
        self.assertClose(result.dev.flux, trueFlux, rtol=0.01)
        self.assertLess(result.dev.getEllipse().getDeterminantRadius(), 0.2)
        self.assertFalse(result.getFlag(result.FAILED))
        self.assertClose(result.flux, trueFlux, rtol=0.01)


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
