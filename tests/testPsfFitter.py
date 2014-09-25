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
import glob
import numpy

import lsst.utils.tests
import lsst.shapelet
import lsst.afw.geom.ellipses
import lsst.meas.multifit
import lsst.meas.base

numpy.random.seed(500)

lsst.pex.logging.Debug("meas.multifit.optimizer.Optimizer", 0)
lsst.pex.logging.Debug("meas.multifit.optimizer.solveTrustRegion", 0)

ELLIPSE_PARAMETER_NAMES = ["eta1", "eta2", "logR", "x", "y"]
DATA_DIR = os.path.join(os.environ["MEAS_MULTIFIT_DIR"], "tests", "data")

def computeMoments(image):
    """Helper function to compute moments of a postage stamp about its origin."""
    maskedImage = lsst.afw.image.MaskedImageD(image)
    result = lsst.meas.base.SdssShapeAlgorithm.Result()
    lsst.meas.base.SdssShapeAlgorithm.apply(
        maskedImage,
        lsst.afw.detection.Footprint(image.getBBox(lsst.afw.image.PARENT)),
        lsst.afw.geom.Point2D(0.0, 0.0),
        result
        )
    return result.getShape()

class PsfFitterTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.configs = {}
        self.configs['fixed'] = lsst.meas.multifit.PsfFitterConfig()
        self.configs['fixed'].primary.ellipticityPriorSigma = 0.0
        self.configs['fixed'].primary.radiusPriorSigma = 0.0
        self.configs['fixed'].primary.positionPriorSigma = 0.0
        self.configs['fixed'].wings.ellipticityPriorSigma = 0.0
        self.configs['fixed'].wings.radiusPriorSigma = 0.0
        self.configs['fixed'].wings.positionPriorSigma = 0.0
        self.configs['ellipse'] = lsst.meas.multifit.PsfFitterConfig()
        self.configs['ellipse'].primary.positionPriorSigma = 0.0
        self.configs['ellipse'].wings.positionPriorSigma = 0.0
        self.configs['full'] = lsst.meas.multifit.PsfFitterConfig()
        self.configs['full'].inner.order = 0
        self.configs['full'].primary.order = 4
        self.configs['full'].wings.order = 4
        self.configs['full'].outer.order = 0
        self.configs['full'].inner.ellipticityPriorSigma = 0.3
        self.configs['full'].inner.radiusPriorSigma = 0.5
        self.configs['full'].inner.positionPriorSigma = 0.1
        self.configs['full'].primary.ellipticityPriorSigma = 0.3
        self.configs['full'].primary.radiusPriorSigma = 0.5
        self.configs['full'].primary.positionPriorSigma = 0.1
        self.configs['full'].wings.ellipticityPriorSigma = 0.3
        self.configs['full'].wings.radiusPriorSigma = 0.5
        self.configs['full'].wings.positionPriorSigma = 0.1
        self.configs['full'].outer.ellipticityPriorSigma = 0.3
        self.configs['full'].outer.radiusPriorSigma = 0.5
        self.configs['full'].outer.positionPriorSigma = 0.1

    def tearDown(self):
        del self.configs

    def testFixedModel(self):
        fitter = lsst.meas.multifit.PsfFitter(self.configs['fixed'].makeControl())
        model = fitter.getModel()

        # check that we have the right numbers and names for parameters
        self.assertEqual(model.getNonlinearDim(), 0)
        self.assertEqual(model.getFixedDim(), 10)
        self.assertEqual(model.getAmplitudeDim(), 2)
        self.assertEqual(model.getBasisCount(), 2)
        self.assertEqual(list(model.getNonlinearNames()), [])
        self.assertEqual(list(model.getAmplitudeNames()), ["primary.alpha[0,0]", "wings.alpha[0,0]"])
        self.assertEqual(list(model.getFixedNames()),
                         ["primary.fiducial.%s" % s for s in ELLIPSE_PARAMETER_NAMES]
                         + ["wings.fiducial.%s" % s for s in ELLIPSE_PARAMETER_NAMES])

        # test that we can round-trip ellipses through the model, and that this agrees
        # with makeShapeletFunction
        ellipseParameters = numpy.array([[0.01, -0.01, 1.1, 0.03, -0.04],
                                         [0.02, -0.02, 0.9, 0.05, -0.06]])
        ellipses1 = model.makeEllipseVector()
        for i in range(len(ellipses1)):
            ellipses1[i].setParameterVector(ellipseParameters[i])
        nonlinear = numpy.zeros(model.getNonlinearDim(), dtype=lsst.meas.multifit.Scalar)
        fixed = numpy.zeros(model.getFixedDim(), dtype=lsst.meas.multifit.Scalar)
        amplitudes = numpy.array([1.0, 0.1], dtype=lsst.meas.multifit.Scalar)
        model.readEllipses(ellipses1, nonlinear, fixed)
        self.assertClose(fixed, ellipseParameters.ravel())
        ellipses2 = model.writeEllipses(nonlinear, fixed)
        msf = model.makeShapeletFunction(nonlinear, amplitudes, fixed)
        self.assertClose(len(msf.getComponents()), len(ellipses1))
        ellipses3 = model.makeEllipseVector()
        for i in range(len(ellipses2)):
            self.assertClose(ellipses1[i].getParameterVector(), ellipses2[i].getParameterVector())
            ellipses3[i] = msf.getComponents()[i].getEllipse()  # need to convert ellipse parametrization
            self.assertClose(ellipses1[i].getParameterVector(), ellipses3[i].getParameterVector())
            self.assertClose(amplitudes[i:i+1], msf.getComponents()[i].getCoefficients())

    def testEllipseModel(self):
        fitter = lsst.meas.multifit.PsfFitter(self.configs['ellipse'].makeControl())
        model = fitter.getModel()

        # check that we have the right numbers and names for parameters
        self.assertEqual(model.getNonlinearDim(), 6)
        self.assertEqual(model.getFixedDim(), 10)
        self.assertEqual(model.getAmplitudeDim(), 2)
        self.assertEqual(model.getBasisCount(), 2)
        self.assertEqual(list(model.getNonlinearNames()),
                         ["primary.%s" % s for s in ELLIPSE_PARAMETER_NAMES[:3]]
                         + ["wings.%s" % s for s in ELLIPSE_PARAMETER_NAMES[:3]]
                         )
        self.assertEqual(list(model.getAmplitudeNames()), ["primary.alpha[0,0]", "wings.alpha[0,0]"])
        self.assertEqual(list(model.getFixedNames()),
                         ["primary.fiducial.%s" % s for s in ELLIPSE_PARAMETER_NAMES]
                         + ["wings.fiducial.%s" % s for s in ELLIPSE_PARAMETER_NAMES])

        # test that we can round-trip ellipses through the model, and that this agrees
        # with makeShapeletFunction
        ellipseParameters = numpy.array([[0.01, -0.01, 1.1, 0.03, -0.04],
                                         [0.02, -0.02, 0.9, 0.05, -0.06]])
        ellipses1 = model.makeEllipseVector()
        for i in range(len(ellipses1)):
            ellipses1[i].setParameterVector(ellipseParameters[i])
        nonlinear = numpy.zeros(model.getNonlinearDim(), dtype=lsst.meas.multifit.Scalar)
        fixed = numpy.zeros(model.getFixedDim(), dtype=lsst.meas.multifit.Scalar)
        amplitudes = numpy.array([1.0, 0.1], dtype=lsst.meas.multifit.Scalar)
        model.readEllipses(ellipses1, nonlinear, fixed)
        self.assertClose(nonlinear, numpy.zeros(model.getNonlinearDim(), dtype=lsst.meas.multifit.Scalar))
        self.assertClose(fixed, ellipseParameters.ravel())
        ellipses2 = model.writeEllipses(nonlinear, fixed)
        msf = model.makeShapeletFunction(nonlinear, amplitudes, fixed)
        self.assertClose(len(msf.getComponents()), len(ellipses1))
        ellipses3 = model.makeEllipseVector()
        for i in range(len(ellipses2)):
            self.assertClose(ellipses1[i].getParameterVector(), ellipses2[i].getParameterVector(), rtol=1E-8)
            ellipses3[i] = msf.getComponents()[i].getEllipse()  # need to convert ellipse parametrization
            self.assertClose(ellipses1[i].getParameterVector(), ellipses3[i].getParameterVector(), rtol=1E-8)
            self.assertClose(amplitudes[i:i+1], msf.getComponents()[i].getCoefficients(), rtol=1E-8)

        # test the ellipse round-tripping again, this time starting with nonzero nonlinear parameters:
        # this will be read back in by adding to the fixed parameters and zeroing the nonlinear parameters.
        nonlinear[:] = 0.5*ellipseParameters[:,:3].ravel()
        ellipses4 = model.writeEllipses(nonlinear, fixed)
        model.readEllipses(ellipses4, nonlinear, fixed)
        self.assertClose(nonlinear, numpy.zeros(model.getNonlinearDim(), dtype=lsst.meas.multifit.Scalar),
                         rtol=1E-8)
        self.assertClose(fixed.reshape(2,5)[:,:3], 1.5*ellipseParameters[:,:3], rtol=1E-8)
        self.assertClose(fixed.reshape(2,5)[:,3:], ellipseParameters[:,3:], rtol=1E-8)

    def testFullModel(self):
        fitter = lsst.meas.multifit.PsfFitter(self.configs['full'].makeControl())
        model = fitter.getModel()

        # check that we have the right numbers and names for parameters
        self.assertEqual(model.getNonlinearDim(), 20)
        self.assertEqual(model.getFixedDim(), 20)
        self.assertEqual(model.getAmplitudeDim(), 2*(1 + lsst.shapelet.computeSize(4)))
        self.assertEqual(model.getBasisCount(), 4)
        self.assertEqual(list(model.getNonlinearNames()),
                         ["inner.%s" % s for s in ELLIPSE_PARAMETER_NAMES]
                         + ["primary.%s" % s for s in ELLIPSE_PARAMETER_NAMES]
                         + ["wings.%s" % s for s in ELLIPSE_PARAMETER_NAMES]
                         + ["outer.%s" % s for s in ELLIPSE_PARAMETER_NAMES]
                         )
        self.assertEqual(list(model.getAmplitudeNames()),
                         ["inner.alpha[0,0]"]
                         + ["primary.alpha[%d,%d]" % (x, y)
                            for n, x, y in lsst.shapelet.HermiteIndexGenerator(4)]
                         + ["wings.alpha[%d,%d]" % (x, y)
                            for n, x, y in lsst.shapelet.HermiteIndexGenerator(4)]
                         + ["outer.alpha[0,0]"])
        self.assertEqual(list(model.getFixedNames()),
                         ["inner.fiducial.%s" % s for s in ELLIPSE_PARAMETER_NAMES]
                         + ["primary.fiducial.%s" % s for s in ELLIPSE_PARAMETER_NAMES]
                         + ["wings.fiducial.%s" % s for s in ELLIPSE_PARAMETER_NAMES]
                         + ["outer.fiducial.%s" % s for s in ELLIPSE_PARAMETER_NAMES]
                         )

        # test that we can round-trip ellipses through the model, and that this agrees
        # with makeShapeletFunction
        ellipseParameters = numpy.array([[0.01, -0.01, 1.1, 0.03, -0.04],
                                         [0.015, -0.015, 1.0, 0.04, -0.05],
                                         [0.02, -0.02, 0.9, 0.05, -0.06],
                                         [0.025, -0.025, 0.8, 0.06, -0.07],
                                         ])
        ellipses1 = model.makeEllipseVector()
        for i in range(len(ellipses1)):
            ellipses1[i].setParameterVector(ellipseParameters[i])
        nonlinear = numpy.zeros(model.getNonlinearDim(), dtype=lsst.meas.multifit.Scalar)
        fixed = numpy.zeros(model.getFixedDim(), dtype=lsst.meas.multifit.Scalar)
        amplitudes = numpy.random.randn(model.getAmplitudeDim())
        model.readEllipses(ellipses1, nonlinear, fixed)
        self.assertClose(nonlinear, numpy.zeros(model.getNonlinearDim(), dtype=lsst.meas.multifit.Scalar))
        self.assertClose(fixed, ellipseParameters.ravel())
        ellipses2 = model.writeEllipses(nonlinear, fixed)
        msf = model.makeShapeletFunction(nonlinear, amplitudes, fixed)
        self.assertClose(len(msf.getComponents()), len(ellipses1))
        ellipses3 = model.makeEllipseVector()
        amplitudeOffset = 0
        for i in range(len(ellipses2)):
            self.assertClose(ellipses1[i].getParameterVector(), ellipses2[i].getParameterVector(), rtol=1E-8)
            ellipses3[i] = msf.getComponents()[i].getEllipse()  # need to convert ellipse parametrization
            amplitudeCount = len(msf.getComponents()[i].getCoefficients())
            self.assertClose(ellipses1[i].getParameterVector(), ellipses3[i].getParameterVector(), rtol=1E-8)
            self.assertClose(amplitudes[amplitudeOffset:amplitudeOffset+amplitudeCount],
                             msf.getComponents()[i].getCoefficients(), rtol=1E-8)
            amplitudeOffset += amplitudeCount

        # test the ellipse round-tripping again, this time starting with nonzero nonlinear parameters:
        # this will be read back in by adding to the fixed parameters and zeroing the nonlinear parameters.
        nonlinear[:] = 0.5*ellipseParameters.ravel()
        ellipses4 = model.writeEllipses(nonlinear, fixed)
        model.readEllipses(ellipses4, nonlinear, fixed)
        self.assertClose(nonlinear, numpy.zeros(model.getNonlinearDim(), dtype=lsst.meas.multifit.Scalar))
        self.assertClose(fixed, 1.5*ellipseParameters.ravel())

    def testApply(self):
        tolerances = {"full": 3E-4, "ellipse": 8E-3, "fixed": 1E-2}
        for filename in glob.glob(os.path.join(DATA_DIR, "psfs", "*.fits")):
            kernelImage = lsst.afw.image.ImageD(filename)
            shape = computeMoments(kernelImage)
            for configKey in ["full", "ellipse", "fixed"]:
                fitter = lsst.meas.multifit.PsfFitter(self.configs[configKey].makeControl())
                multiShapeletFit = fitter.apply(kernelImage, shape, 0.01)
                modelImage = lsst.afw.image.ImageD(kernelImage.getBBox(lsst.afw.image.PARENT))
                multiShapeletFit.evaluate().addToImage(modelImage)
                self.assertClose(kernelImage.getArray(), modelImage.getArray(),
                                 atol=tolerances[configKey],
                                 plotOnFailure=True)

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(PsfFitterTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
