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
import numpy

import lsst.pex.logging
import lsst.utils.tests
import lsst.shapelet.tests
import lsst.afw.geom.ellipses
import lsst.afw.image
import lsst.afw.math
import lsst.afw.detection
import lsst.meas.modelfit
import lsst.meas.modelfit.display
import lsst.afw.display.ds9

numpy.random.seed(500)

ASSERT_CLOSE_KWDS = dict(plotOnFailure=True, printFailures=False)

def makeGaussianFunction(ellipse, flux=1.0):
    """Create a single-Gaussian MultiShapeletFunction

    ellipse may be an afw.geom.ellipses.Ellipse or a float radius for a circle
    """
    s = lsst.shapelet.ShapeletFunction(0, lsst.shapelet.HERMITE, ellipse)
    s.getCoefficients()[0] = 1.0
    s.normalize()
    s.getCoefficients()[0] *= flux
    msf = lsst.shapelet.MultiShapeletFunction()
    msf.getComponents().push_back(s)
    return msf

def addGaussian(exposure, ellipse, flux, psf=None):
    s = makeGaussianFunction(ellipse, flux)
    if psf is not None:
        s = s.convolve(psf)
    imageF = exposure.getMaskedImage().getImage()
    imageD = lsst.afw.image.ImageD(imageF.getBBox())
    s.evaluate().addToImage(imageD)
    imageF += imageD.convertF()

def scaleExposure(exposure, factor):
    mi = exposure.getMaskedImage()
    mi *= factor

class UnitTransformedLikelihoodTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.position = lsst.afw.coord.IcrsCoord(45.0*lsst.afw.geom.degrees, 45.0*lsst.afw.geom.degrees)
        self.model = lsst.meas.modelfit.Model.makeGaussian(lsst.meas.modelfit.Model.FIXED_CENTER)
        self.ellipse = lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes(6.0, 5.0, numpy.pi/6))
        self.flux = 50.0
        ev = self.model.makeEllipseVector()
        ev[0] = self.ellipse
        self.nonlinear = numpy.zeros(self.model.getNonlinearDim(), dtype=lsst.meas.modelfit.Scalar)
        self.fixed = numpy.zeros(self.model.getFixedDim(), dtype=lsst.meas.modelfit.Scalar)
        self.model.readEllipses(ev, self.nonlinear, self.fixed)
        self.amplitudes = numpy.zeros(self.model.getAmplitudeDim(), dtype=lsst.meas.modelfit.Scalar)
        self.amplitudes[:] = self.flux
        # setup ideal exposure0: uses fit Wcs and Calib, has delta function PSF
        cdelt = (0.2*lsst.afw.geom.arcseconds).asDegrees()
        wcs0 = lsst.afw.image.makeWcs(self.position, lsst.afw.geom.Point2D(), cdelt, 0.0, 0.0, cdelt)
        calib0 = lsst.afw.image.Calib()
        calib0.setFluxMag0(10000)
        self.psf0 = makeGaussianFunction(0.0)
        self.bbox0 = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-100, -100), lsst.afw.geom.Point2I(100, 100))
        self.footprint0 = lsst.afw.detection.Footprint(self.bbox0)
        self.exposure0 = lsst.afw.image.ExposureF(self.bbox0)
        self.exposure0.setWcs(wcs0)
        self.exposure0.setCalib(calib0)
        self.sys0 = lsst.meas.modelfit.UnitSystem(self.exposure0)
        addGaussian(self.exposure0, self.ellipse, self.flux, psf=self.psf0)
        self.exposure0.getMaskedImage().getVariance().set(1.0)
        # setup secondary exposure: 2x pixel scale, 3x gain, Gaussian PSF with sigma=2.5pix
        cdelt = (0.4*lsst.afw.geom.arcseconds).asDegrees()
        wcs1 = lsst.afw.image.makeWcs(self.position, lsst.afw.geom.Point2D(), cdelt, 0.0, 0.0, cdelt)
        calib1 = lsst.afw.image.Calib()
        calib1.setFluxMag0(30000)
        self.sys1 = lsst.meas.modelfit.UnitSystem(wcs1, calib1)
        # transform object that maps between exposures (not including PSF)
        self.t01 = lsst.meas.modelfit.LocalUnitTransform(self.position, self.sys0, self.sys1)
        self.bbox1 = lsst.afw.geom.Box2I(self.bbox0)
        self.bbox1.grow(-60)
        self.footprint1 = lsst.afw.detection.Footprint(self.bbox1)
        self.psfSigma1 = 2.5
        self.psf1 = makeGaussianFunction(self.psfSigma1)

    def tearDown(self):
        del self.position
        del self.model
        del self.ellipse
        del self.bbox0
        del self.footprint0
        del self.exposure0
        del self.sys0
        del self.sys1
        del self.t01
        del self.bbox1
        del self.footprint1
        del self.psf1

    def checkLikelihood(self, likelihood, data):
        self.assertClose(likelihood.getData().reshape(data.shape), data, rtol=1E-6, **ASSERT_CLOSE_KWDS)
        matrix = numpy.zeros((1, likelihood.getDataDim()), dtype=lsst.meas.modelfit.Pixel).transpose()
        likelihood.computeModelMatrix(matrix, self.nonlinear)
        model = numpy.dot(matrix, self.amplitudes)
        self.assertClose(model.reshape(data.shape), data, rtol=1E-6, atol=1E-7, **ASSERT_CLOSE_KWDS)

    def testModel(self):
        """Test that when we use a Model to create a MultiShapeletFunction from our parameter vectors
        it agrees with the reimplementation here."""
        msf = self.model.makeShapeletFunction(self.nonlinear, self.amplitudes, self.fixed)
        image0a = lsst.afw.image.ImageD(self.bbox0)
        msf.evaluate().addToImage(image0a)
        self.assertClose(image0a.getArray(), self.exposure0.getMaskedImage().getImage().getArray(),
                         rtol=1E-6, atol=1E-7, **ASSERT_CLOSE_KWDS)

    def testWarp(self):
        """Test that transforming ellipses and fluxes with LocalUnitTransform agrees with warping
        """
        warpCtrl = lsst.afw.math.WarpingControl("lanczos5")
        # exposure1: check image; just a transform and scaling of exposure0 (no PSF...yet)
        exposure1 = lsst.afw.image.ExposureF(self.bbox1)
        addGaussian(exposure1, self.ellipse.transform(self.t01.geometric), self.flux * self.t01.flux)
        # exposure1a: warp exposure0 using warpExposure with WCS arguments
        exposure1a = lsst.afw.image.ExposureF(self.bbox1)
        exposure1a.setWcs(self.sys1.wcs)
        lsst.afw.math.warpExposure(exposure1a, self.exposure0, warpCtrl)
        exposure1a.setCalib(self.sys1.calib)
        scaleExposure(exposure1a, self.t01.flux)
        self.assertClose(exposure1.getMaskedImage().getImage().getArray(),
                         exposure1a.getMaskedImage().getImage().getArray(),
                         rtol=1E-6, **ASSERT_CLOSE_KWDS)
        # exposure1b: warp exposure0 using warpImage with AffineTransform arguments
        exposure1b = lsst.afw.image.ExposureF(self.bbox1)
        exposure1b.setWcs(self.sys1.wcs)
        lsst.afw.math.warpImage(exposure1b.getMaskedImage(), self.exposure0.getMaskedImage(),
                                self.t01.geometric, warpCtrl)
        exposure1b.setCalib(self.sys1.calib)
        scaleExposure(exposure1b, self.t01.flux)
        self.assertClose(exposure1.getMaskedImage().getImage().getArray(),
                         exposure1b.getMaskedImage().getImage().getArray(),
                         rtol=1E-6, **ASSERT_CLOSE_KWDS)
        # now we rebuild exposure1 with the PSF convolution included, and convolve 1a->1c using an
        # afw::math::Kernel.  Since 1a is the same as 1b, there's no need to convolve 1b too.
        exposure1 = lsst.afw.image.ExposureF(self.bbox1)
        addGaussian(exposure1, self.ellipse.transform(self.t01.geometric), self.flux * self.t01.flux,
                    psf=self.psf1)
        kernel = lsst.afw.math.AnalyticKernel(
            int(self.psfSigma1*16)+1, int(self.psfSigma1*16)+1,
            lsst.afw.math.GaussianFunction2D(self.psfSigma1, self.psfSigma1)
            )
        exposure1c = lsst.afw.image.ExposureF(self.bbox1)
        ctrl = lsst.afw.math.ConvolutionControl()
        ctrl.setDoCopyEdge(True)
        lsst.afw.math.convolve(exposure1c.getMaskedImage(), exposure1a.getMaskedImage(), kernel, ctrl)
        self.assertClose(exposure1.getMaskedImage().getImage().getArray(),
                         exposure1c.getMaskedImage().getImage().getArray(),
                         rtol=1E-5, atol=1E-6, **ASSERT_CLOSE_KWDS)

    def testDirect(self):
        """Test likelihood evaluation when the fit system is the same as the data system.
        """
        ctrl = lsst.meas.modelfit.UnitTransformedLikelihoodControl()
        var = numpy.random.rand(self.bbox0.getHeight(), self.bbox0.getWidth()) + 2.0
        self.exposure0.getMaskedImage().getVariance().getArray()[:,:] = var
        efv = lsst.meas.modelfit.EpochFootprintVector()
        efv.push_back(lsst.meas.modelfit.EpochFootprint(self.footprint0, self.exposure0, self.psf0))
        # test with per-pixel weights, using both ctors
        ctrl.usePixelWeights = True
        data = self.exposure0.getMaskedImage().getImage().getArray() / var**0.5
        l0a = lsst.meas.modelfit.UnitTransformedLikelihood(self.model, self.fixed, self.sys0, self.position,
                                                     self.exposure0, self.footprint0, self.psf0, ctrl)
        self.checkLikelihood(l0a, data)
        l0b = lsst.meas.modelfit.UnitTransformedLikelihood(self.model, self.fixed, self.sys0, self.position,
                                                     efv, ctrl)
        self.checkLikelihood(l0b, data)
        # test with constant weights, using both ctors
        ctrl.usePixelWeights = False
        data = self.exposure0.getMaskedImage().getImage().getArray()
        l0c = lsst.meas.modelfit.UnitTransformedLikelihood(self.model, self.fixed, self.sys0, self.position,
                                                           self.exposure0, self.footprint0, self.psf0, ctrl)
        self.checkLikelihood(l0c, data)
        l0d = lsst.meas.modelfit.UnitTransformedLikelihood(self.model, self.fixed, self.sys0, self.position,
                                                     efv, ctrl)
        self.checkLikelihood(l0d, data)

    def testProjected(self):
        """Test likelihood evaluation when the fit system is not the same as the data system.
        """
        # Start by building the data exposure
        exposure1 = lsst.afw.image.ExposureF(self.bbox1)
        addGaussian(exposure1, self.ellipse.transform(self.t01.geometric), self.flux * self.t01.flux,
                    psf=self.psf1)
        exposure1.setWcs(self.sys1.wcs)
        exposure1.setCalib(self.sys1.calib)
        var = numpy.random.rand(self.bbox1.getHeight(), self.bbox1.getWidth()) + 2.0
        exposure1.getMaskedImage().getVariance().getArray()[:,:] = var
        ctrl = lsst.meas.modelfit.UnitTransformedLikelihoodControl()
        efv = lsst.meas.modelfit.EpochFootprintVector()
        efv.push_back(lsst.meas.modelfit.EpochFootprint(self.footprint1, exposure1, self.psf1))
        # test with per-pixel weights, using both ctors
        ctrl.usePixelWeights = True
        data = exposure1.getMaskedImage().getImage().getArray() / var**0.5
        l1a = lsst.meas.modelfit.UnitTransformedLikelihood(self.model, self.fixed, self.sys0, self.position,
                                                     exposure1, self.footprint1, self.psf1, ctrl)
        self.checkLikelihood(l1a, data)
        l1b = lsst.meas.modelfit.UnitTransformedLikelihood(self.model, self.fixed, self.sys0, self.position,
                                                     efv, ctrl)
        self.checkLikelihood(l1b, data)
        # test with constant weights, using both ctors
        ctrl.usePixelWeights = False
        data = exposure1.getMaskedImage().getImage().getArray()
        l1c = lsst.meas.modelfit.UnitTransformedLikelihood(self.model, self.fixed, self.sys0, self.position,
                                                     exposure1, self.footprint1, self.psf1, ctrl)
        self.checkLikelihood(l1c, data)
        l1d = lsst.meas.modelfit.UnitTransformedLikelihood(self.model, self.fixed, self.sys0, self.position,
                                                     efv, ctrl)
        self.checkLikelihood(l1d, data)

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(UnitTransformedLikelihoodTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
