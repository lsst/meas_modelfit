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

log = lsst.pex.logging.Debug("meas.multifit.optimizer.Optimizer", 0)
log = lsst.pex.logging.Debug("meas.multifit.optimizer.solveTrustRegion", 0)

class PsfFittingTestCase(lsst.utils.tests.TestCase):

    def testMultiShapeletPsfModel(self):
        orders = [2, 3, 4]
        model = lsst.meas.multifit.makeMultiShapeletPsfModel(orders)
        self.assertEqual(model.getNonlinearDim(), 5*len(orders))
        self.assertEqual(model.getAmplitudeDim(), sum(lsst.shapelet.computeSize(order) for order in orders))
        self.assertEqual(model.getFixedDim(), 0)
        self.assertEqual(model.getBasisCount(), len(orders))
        nonlinearNames = ["%d.%s" % (n, s) for n in range(len(orders)) for s in ["eta1", "eta2", "logR"]]
        nonlinearNames += ["%d.%s" % (n, s) for n in range(len(orders)) for s in ["x", "y"]]
        self.assertEqual(list(model.getNonlinearNames()), nonlinearNames)
        amplitudeNames = []
        for i, order in enumerate(orders):
            amplitudeNames.extend("%d.alpha%d" % (i, n) for n in range(lsst.shapelet.computeSize(order)))
        self.assertEqual(list(model.getAmplitudeNames()), amplitudeNames)
        self.assertEqual(list(model.getFixedNames()), [])

    def testSingleGaussianFit(self):
        orders = [0]
        sigma = 0.001
        model = lsst.meas.multifit.makeMultiShapeletPsfModel(orders)
        nonlinear = numpy.zeros(model.getNonlinearDim(), dtype=lsst.meas.multifit.Scalar)
        amplitudes = numpy.ones(model.getAmplitudeDim(), dtype=lsst.meas.multifit.Scalar)
        fixed = numpy.zeros(model.getFixedDim(), dtype=lsst.meas.multifit.Scalar)
        ellipses = model.makeEllipseVector()
        ellipses[0].setCore(lsst.afw.geom.ellipses.Axes(4.0, 4.0, 0.0))
        model.readEllipses(ellipses, nonlinear, fixed)
        msf = model.makeShapeletFunction(nonlinear, amplitudes, fixed)
        image = numpy.zeros((25, 25), dtype=float)
        xy0 = lsst.afw.geom.Point2I(-12, -12)
        msf.evaluate().addToImage(image, xy0)
        image = image.astype(lsst.meas.multifit.Pixel)
        likelihood = lsst.meas.multifit.MultiShapeletPsfLikelihood(image, xy0, model, sigma, fixed)
        self.assertEqual(likelihood.getDataDim(), image.size)
        self.assertEqual(likelihood.getAmplitudeDim(), model.getAmplitudeDim())
        self.assertEqual(likelihood.getNonlinearDim(), model.getNonlinearDim())
        self.assertEqual(likelihood.getFixedDim(), model.getFixedDim())
        self.assertClose(likelihood.getData()*sigma, image.flat, rtol=1E-6)
        matrix = numpy.zeros((likelihood.getAmplitudeDim(), likelihood.getDataDim()),
                             dtype=lsst.meas.multifit.Pixel).transpose()
        likelihood.computeModelMatrix(matrix, nonlinear)
        self.assertClose(image, matrix.reshape(image.shape)*sigma, rtol=5E-7, atol=1E-7, plotOnFailure=True)
        parameters = numpy.concatenate([nonlinear, amplitudes])
        objective = lsst.meas.multifit.OptimizerObjective.makeFromLikelihood(likelihood, None)
        ctrl = lsst.meas.multifit.OptimizerControl()
        optimizer = lsst.meas.multifit.Optimizer(objective, parameters, ctrl)
        n = optimizer.run()
        self.assert_(optimizer.getState() & optimizer.CONVERGED)
        self.assertClose(optimizer.getParameters(), parameters, atol=1E-6, rtol=1E-8)
        print "iterations=%d, state=0x%05x, parameters=%s" % (n, optimizer.getState(),
                                                              optimizer.getParameters())
        for i in range(10):
            parameters0 = parameters + numpy.random.randn(parameters.size) * 0.1
            optimizer = lsst.meas.multifit.Optimizer(objective, parameters0, ctrl)
            n = optimizer.run()
            self.assert_(optimizer.getState() & optimizer.CONVERGED)
            self.assertClose(optimizer.getParameters(), parameters, atol=5E-7, rtol=1E-8)
            print "iterations=%d, state=0x%05x, parameters=%s" % (n, optimizer.getState(),
                                                                  optimizer.getParameters())

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(PsfFittingTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
