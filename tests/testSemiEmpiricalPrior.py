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

lsst.pex.logging.Debug("meas.modelfit.SemiEmpiricalPrior", 10)


class SemiEmpiricalPriorTestCase(lsst.utils.tests.TestCase):

    NUM_DIFF_STEP = 1E-4

    def setUp(self):
        # a prior with broad ramps and non-zero slope; broad ramps makes evaluating numerical
        # derivatives easier, and we want to do that to check the analytic ones
        self.ctrl = lsst.meas.modelfit.SemiEmpiricalPrior.Control()
        self.ctrl.ellipticityCore = 4.0
        self.ctrl.ellipticitySigma = 10.0
        self.ctrl.logRadiusMinOuter = self.ctrl.logRadiusMinInner - 8.0
        self.ctrl.logRadiusMu = 2.0
        self.ctrl.logRadiusSigma = 5.0
        self.ctrl.logRadiusNu = 2.0
        self.prior = lsst.meas.modelfit.SemiEmpiricalPrior(self.ctrl)
        self.amplitudes = numpy.array([1.0], dtype=lsst.meas.modelfit.Scalar)
        dtype = numpy.dtype([("eta1", float), ("eta2", float), ("lnR", float), ("p", float),
                             ("d_eta1", float), ("d_eta2", float), ("d_lnR", float),
                             ("d2_eta1_eta1", float), ("d2_eta1_eta2", float),
                             ("d2_eta1_lnR", float), ("d2_eta2_eta2", float),
                             ("d2_eta2_lnR", float), ("d2_lnR_lnR", float)])
        self.data = numpy.loadtxt("tests/data/SEP.txt", dtype=dtype)

    def tearDown(self):
        del self.prior
        del self.amplitudes

    def testEvaluate(self):
        for row in self.data:
            p = self.prior.evaluate(numpy.array([row["eta1"], row["eta2"], row["lnR"]]), self.amplitudes)
            self.assertClose(row["p"], p)

    def testGradient(self):
        for row in self.data:
            grad = numpy.zeros(4, dtype=float)
            hess = numpy.zeros((4, 4), dtype=float)
            self.prior.evaluateDerivatives(
                numpy.array([row["eta1"], row["eta2"], row["lnR"]]),
                self.amplitudes,
                grad[:3], grad[3:],
                hess[:3, :3], hess[3:, 3:], hess[:3, 3:]
            )
            self.assertClose(row["d_eta1"], grad[0])
            self.assertClose(row["d_eta2"], grad[1])
            self.assertClose(row["d_lnR"], grad[2])

    def testHessian(self):
        for row in self.data:
            grad = numpy.zeros(4, dtype=float)
            hess = numpy.zeros((4, 4), dtype=float)
            self.prior.evaluateDerivatives(
                numpy.array([row["eta1"], row["eta2"], row["lnR"]]),
                self.amplitudes,
                grad[:3], grad[3:],
                hess[:3, :3], hess[3:, 3:], hess[:3, 3:]
            )
            self.assertClose(row["d2_eta1_eta1"], hess[0, 0])
            self.assertClose(row["d2_eta1_eta2"], hess[0, 1])
            self.assertClose(row["d2_eta1_lnR"], hess[0, 2])
            self.assertClose(row["d2_eta2_eta2"], hess[1, 1])
            self.assertClose(row["d2_eta2_lnR"], hess[1, 2])
            self.assertClose(row["d2_lnR_lnR"], hess[2, 2])

    def evaluatePrior(self, eta1, eta2, lnR):
        b = numpy.broadcast(eta1, eta2, lnR)
        p = numpy.zeros(b.shape, dtype=lsst.meas.modelfit.Scalar)
        for i, (eta1i, eta2i, lnRi) in enumerate(b):
            p.flat[i] = self.prior.evaluate(numpy.array([eta1i, eta2i, lnRi]), self.amplitudes)
        return p


def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(SemiEmpiricalPriorTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)


def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
