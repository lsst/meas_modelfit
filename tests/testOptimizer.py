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
import lsst.pex.logging
import lsst.meas.multifit

numpy.random.seed(500)

log = lsst.pex.logging.Debug("meas.multifit.optimizer", 10)

class OptimizerTestCase(lsst.utils.tests.TestCase):

    def testTrustRegionSolver(self):
        tolerance = 1E-6
        # start with some positive definite matrices, constructed from random least-squares problems
        log.info("Testing solveTrustRegion with positive-definite matrices")
        m = numpy.random.randn(30, 5)
        y = numpy.random.randn(30)
        f = numpy.dot(m.transpose(), m)
        g = numpy.dot(m.transpose(), y)
        x = numpy.zeros(5)
        for r in numpy.linspace(1E-3, 0.8, 5):
            lsst.meas.multifit.solveTrustRegion(x, f, g, r, tolerance)
            self.assertLessEqual(numpy.linalg.norm(x), r * (1.0 + tolerance))
        # now we try some matrices with zero eigenvalues due to model degeneracies
        log.info("Testing solveTrustRegion with positive-semidefinite matrices")
        m[:,-1] = m[:,0]
        f = numpy.dot(m.transpose(), m)
        g = numpy.dot(m.transpose(), y)
        for r in numpy.linspace(1E-3, 0.8, 5):
            lsst.meas.multifit.solveTrustRegion(x, f, g, r, tolerance)
            self.assertLessEqual(numpy.linalg.norm(x), r * (1.0 + tolerance))
        m[:,-2] = m[:,1]
        f = numpy.dot(m.transpose(), m)
        g = numpy.dot(m.transpose(), y)
        for r in numpy.linspace(1E-3, 0.8, 5):
            lsst.meas.multifit.solveTrustRegion(x, f, g, r, tolerance)
            self.assertLessEqual(numpy.linalg.norm(x), r * (1.0 + tolerance))
        log.info("Testing solveTrustRegion with indefinite matrices")
        for i in range(3):
            m = numpy.random.randn(5, 5)
            f = m + m.transpose()
            g = numpy.random.randn(5)
            for r in numpy.linspace(1E-3, 0.8, 5):
                lsst.meas.multifit.solveTrustRegion(x, f, g, r, tolerance)
                self.assertLessEqual(numpy.linalg.norm(x), r * (1.0 + tolerance))

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(OptimizerTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
