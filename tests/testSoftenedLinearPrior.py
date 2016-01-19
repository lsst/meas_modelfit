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
import lsst.meas.modelfit

try:
    import scipy.integrate
except ImportError:
    scipy = None

numpy.random.seed(500)

lsst.pex.logging.Debug("meas.modelfit.SoftenedLinearPrior", 10)

class SoftenedLinearPriorTestCase(lsst.utils.tests.TestCase):

    NUM_DIFF_STEP = 1E-3

    def setUp(self):
        # a prior with broad ramps and non-zero slope; broad ramps makes evaluating numerical
        # derivatives easier, and we want to do that to check the analytic ones
        ctrl = lsst.meas.modelfit.SoftenedLinearPrior.Control()
        ctrl.logRadiusMinOuter = ctrl.logRadiusMinInner - 2.0
        ctrl.logRadiusMaxOuter = ctrl.logRadiusMaxInner + 2.0
        ctrl.ellipticityMaxOuter = ctrl.ellipticityMaxInner + 2.0
        ctrl.logRadiusMinMaxRatio = 2.0
        self.prior = lsst.meas.modelfit.SoftenedLinearPrior(ctrl)
        self.amplitudes = numpy.array([1.0], dtype=lsst.meas.modelfit.Scalar)

    def tearDown(self):
        del self.prior
        del self.amplitudes

    def evaluatePrior(self, e1, e2, r):
        b = numpy.broadcast(e1, e2, r)
        p = numpy.zeros(b.shape, dtype=lsst.meas.modelfit.Scalar)
        for i, (e1i, e2i, ri) in enumerate(b):
            p.flat[i] = self.prior.evaluate(numpy.array([e1i, e2i, ri]), self.amplitudes)
        return p

    def checkDerivatives(self, e1, e2, r):
        nonlinear = numpy.array([e1, e2, r], dtype=lsst.meas.modelfit.Scalar)
        amplitudeGradient = numpy.zeros(1, dtype=lsst.meas.modelfit.Scalar)
        amplitudeHessian = numpy.zeros((1,1), dtype=lsst.meas.modelfit.Scalar)
        crossHessian = numpy.zeros((3,1), dtype=lsst.meas.modelfit.Scalar)
        nonlinearGradient = numpy.zeros(3, dtype=lsst.meas.modelfit.Scalar)
        nonlinearHessian = numpy.zeros((3, 3), dtype=lsst.meas.modelfit.Scalar)
        self.prior.evaluateDerivatives(nonlinear, self.amplitudes,
                                  nonlinearGradient, amplitudeGradient,
                                  nonlinearHessian, amplitudeHessian,
                                  crossHessian)
        p = self.prior.evaluate(nonlinear, self.amplitudes)
        for i in range(3):
            nonlinearA = nonlinear.copy()
            nonlinearB = nonlinear.copy()
            nonlinearA[i] -= self.NUM_DIFF_STEP
            nonlinearB[i] += self.NUM_DIFF_STEP
            pA = self.prior.evaluate(nonlinearA, self.amplitudes)
            pB = self.prior.evaluate(nonlinearB, self.amplitudes)
            dp = (pB - pA) / (2*self.NUM_DIFF_STEP)
            self.assertClose(nonlinearGradient[i], dp, rtol=1E-3, atol=1E-8)
            d2p = (pA + pB - 2*p) / self.NUM_DIFF_STEP**2
            self.assertClose(nonlinearHessian[i,i], d2p, rtol=1E-3, atol=1E-8)
            for j in range(i+1, 3):
                nonlinearAA = nonlinearA.copy()
                nonlinearAB = nonlinearA.copy()
                nonlinearBA = nonlinearB.copy()
                nonlinearBB = nonlinearB.copy()
                nonlinearAA[j] -= self.NUM_DIFF_STEP
                nonlinearAB[j] += self.NUM_DIFF_STEP
                nonlinearBA[j] -= self.NUM_DIFF_STEP
                nonlinearBB[j] += self.NUM_DIFF_STEP
                pAA = self.prior.evaluate(nonlinearAA, self.amplitudes)
                pAB = self.prior.evaluate(nonlinearAB, self.amplitudes)
                pBA = self.prior.evaluate(nonlinearBA, self.amplitudes)
                pBB = self.prior.evaluate(nonlinearBB, self.amplitudes)
                d2p = (pBB - pAB - pBA + pAA) / (4*self.NUM_DIFF_STEP**2)
                self.assertClose(nonlinearHessian[i,j], d2p, rtol=1E-3, atol=1E-8)

    def testDerivatives(self):
        """Test that evaluateDerivatives() returns results similar to finite-differences
        on evaluate().
        """
        ctrl = self.prior.getControl()
        # a single |e| value for each ellipticity zone
        ellipticityPoints = numpy.array([0.5*ctrl.ellipticityMaxInner,
                                         0.5*(ctrl.ellipticityMaxInner + ctrl.ellipticityMaxOuter)])
        # a single ln(radius) value for each logRadius zone
        logRadiusPoints = numpy.array([0.5*(ctrl.logRadiusMinOuter + ctrl.logRadiusMinInner),
                                       0.5*(ctrl.logRadiusMinInner + ctrl.logRadiusMaxInner),
                                       0.5*(ctrl.logRadiusMaxInner + ctrl.logRadiusMaxOuter)])
        # a range of position angles
        thetaPoints = numpy.linspace(0.0, numpy.pi, 5)
        for theta in thetaPoints:
            for e in ellipticityPoints:
                e1 = e*numpy.cos(2.0*theta)
                e2 = e*numpy.sin(2.0*theta)
                for r in logRadiusPoints:
                    self.checkDerivatives(e1, e2, r)

    @unittest.skipIf(scipy is None, "could not import scipy")
    def testIntegral(self):
        """Test that the prior is properly normalized.

        Normally, this test has a very low bar, because it's expensive to compute a high-quality
        numerical integral to compare with.  Even so, the scipy integrator does better than it
        thinks it does, and we use that smaller tolerance for the test.  That means this test
        could fail if something about the scipy integrator changes, because we're not telling it
        that it has to get as close as it currently is (because doing so would take way too long).

        If this class is ever changed, we should do at least one of this test with the tolerances
        tightened.
        """
        ctrl = self.prior.getControl()
        integral, absErr = scipy.integrate.tplquad(
            self.evaluatePrior,
            ctrl.logRadiusMinOuter, ctrl.logRadiusMaxOuter,
            lambda logR: -ctrl.ellipticityMaxOuter,
            lambda logR: ctrl.ellipticityMaxOuter,
            lambda logR, e2: -(ctrl.ellipticityMaxOuter**2 - e2**2)**0.5,
            lambda logR, e2: (ctrl.ellipticityMaxOuter**2 - e2**2)**0.5,
            epsabs=1.0,
            epsrel=1.0,
            )
        self.assertClose(integral, 1.0, atol=0.01)

    def testEllipticityDistribution(self):
        """Test that the ellipticity distribution is constant in the inner region,
        mononotically decreasing in the ramp, and zero in the outer region, according
        to evaluate().
        """
        ctrl = self.prior.getControl()
        # a range of |e| values in each ellipticity zone
        eInnerPoints = numpy.linspace(0.0, ctrl.ellipticityMaxInner, 5)
        eRampPoints = numpy.linspace(ctrl.ellipticityMaxInner, ctrl.ellipticityMaxOuter, 5)
        eOuterPoints = numpy.linspace(ctrl.ellipticityMaxOuter, ctrl.ellipticityMaxOuter + 5.0, 5)
        # a range of position angles
        thetaPoints = numpy.linspace(0.0, numpy.pi, 5)
        # a single ln(radius) value for each logRadius zone
        logRadiusPoints = numpy.array([0.5*(ctrl.logRadiusMinOuter + ctrl.logRadiusMinInner),
                                       0.5*(ctrl.logRadiusMinInner + ctrl.logRadiusMaxInner),
                                       0.5*(ctrl.logRadiusMaxInner + ctrl.logRadiusMaxOuter)])
        for logRadius in logRadiusPoints:
            for theta in thetaPoints:
                # All inner points should have the same value
                pInner = self.evaluatePrior(eInnerPoints*numpy.cos(2*theta),
                                            eInnerPoints*numpy.sin(2*theta),
                                            logRadius)
                self.assertClose(pInner.mean(), pInner)

                # Each ramp point should be greater than the next one
                pRamp = self.evaluatePrior(eRampPoints*numpy.cos(2*theta),
                                           eRampPoints*numpy.sin(2*theta),
                                           logRadius)
                self.assertTrue((pRamp[:-1] > pRamp[1:]).all())

                # Each outer point should be zero
                pOuter = self.evaluatePrior(eOuterPoints*numpy.cos(2*theta),
                                            eOuterPoints*numpy.sin(2*theta),
                                            logRadius)
                self.assertClose(pOuter, 0.0)

    def testLogRadiusDistribution(self):
        """Test that the ln(radius) distribution is constant in the inner region,
        mononotically decreasing in the ramps, and zero in the outer regions, according
        to evaluate().
        """
        ctrl = self.prior.getControl()
        # a range of ln(radius) values in each logRadius zone
        rLowerOuterPoints = numpy.linspace(ctrl.logRadiusMinOuter - 2.0, ctrl.logRadiusMinOuter, 5)
        rLowerRampPoints = numpy.linspace(ctrl.logRadiusMinOuter, ctrl.logRadiusMinInner, 5)
        rInnerPoints = numpy.linspace(ctrl.logRadiusMinInner, ctrl.logRadiusMaxInner, 5)
        rUpperRampPoints = numpy.linspace(ctrl.logRadiusMaxInner, ctrl.logRadiusMaxOuter, 5)
        rUpperOuterPoints = numpy.linspace(ctrl.logRadiusMaxOuter, ctrl.logRadiusMaxOuter + 2.0, 5)

        # a range of position angles
        thetaPoints = numpy.linspace(0.0, numpy.pi, 5)
        # a single |e| value for each ellipticity zone
        ellipticityPoints = numpy.array([0.5*ctrl.ellipticityMaxInner,
                                         0.5*(ctrl.ellipticityMaxInner + ctrl.ellipticityMaxOuter)])
        for ellipticity in ellipticityPoints:
            for theta in thetaPoints:
                e1 = ellipticity*numpy.cos(2*theta)
                e2 = ellipticity*numpy.sin(2*theta)
                # Outer points should be zero
                pLowerOuter = self.evaluatePrior(e1, e2, rLowerOuterPoints)
                self.assertClose(pLowerOuter, 0.0)
                # Each ramp point should be less than the next one
                pLowerRamp = self.evaluatePrior(e1, e2, rLowerRampPoints)
                self.assertTrue((pLowerRamp[:-1] < pLowerRamp[1:]).all())
                # All adjacent inner points should have the same distance between them (constant slope)
                pInner = self.evaluatePrior(e1, e2, rInnerPoints)
                diffs = pInner[1:] - pInner[:-1]
                self.assertClose(diffs.mean(), diffs)
                # Each ramp point should be greater than the next one
                pUpperRamp = self.evaluatePrior(e1, e2, rUpperRampPoints)
                self.assertTrue((pUpperRamp[:-1] > pUpperRamp[1:]).all())
                # Outer points should be zero
                pUpperOuter = self.evaluatePrior(e1, e2, rUpperOuterPoints)
                self.assertClose(pUpperOuter, 0.0)


def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(SoftenedLinearPriorTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
