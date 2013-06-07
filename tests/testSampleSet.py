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

import lsst.utils.tests
import lsst.afw.geom.ellipses
import lsst.afw.image
import lsst.afw.detection
import lsst.shapelet.tests
import lsst.meas.multifit

numpy.random.seed(500)

class SampleSetTestCase(lsst.shapelet.tests.ShapeletTestCase):

    def logDist(self, p):
        return 0.5*numpy.sum(((p - self.mu) / self.sigma)**2) + self.norm1

    def setUp(self):
        point = lsst.meas.multifit.SamplePoint(2, 1)
        point.joint.mu = numpy.array([25.0], dtype=numpy.float32)
        point.joint.fisher = numpy.array([[1.0]], dtype=numpy.float32)
        self.samples = lsst.meas.multifit.SampleSet(2, 1, "SeparableConformalShearLogTraceRadius")
        self.mu = numpy.array([0.15, -0.12], dtype=float)
        self.sigma = numpy.array([1.1, 0.9], dtype=float)
        self.norm1 = -numpy.log(numpy.product(self.sigma)/(2.0*numpy.pi))
        self.norm2 = 0.5*numpy.log(2.0*numpy.pi)
        for n in range(20000):
            point.parameters = (numpy.random.randn(2) * self.sigma + self.mu).astype(numpy.float32)
            point.proposal = self.logDist(point.parameters)
            point.joint.r = point.proposal + self.norm2
            self.samples.add(point)

    def tearDown(self):
        del self.samples
        del self.mu
        del self.sigma

    def testEstimators(self):
        logSum1 = self.samples.applyPrior(lsst.meas.multifit.FlatPrior())
        cat = self.samples.asCatalog()
        self.assertClose(cat['marginal'], cat['proposal'], rtol=1E-6)
        self.assertClose(cat['weight'], 1.0 / self.samples.size(), rtol=1E-6)
        self.assertClose(numpy.exp(-logSum1), 1.0, rtol=1E-7)
        self.assertClose(self.samples.computeMean(), numpy.mean(cat['parameters'], axis=0), rtol=1E-5)
        quantiles = self.samples.computeQuantiles(numpy.array([0.2, 0.4, 0.6, 0.8]))
        self.assertClose(quantiles[0,:], numpy.percentile(cat['parameters'], 20, axis=0), rtol=2E-3)
        self.assertClose(quantiles[1,:], numpy.percentile(cat['parameters'], 40, axis=0), rtol=2E-3)
        self.assertClose(quantiles[2,:], numpy.percentile(cat['parameters'], 60, axis=0), rtol=2E-3)
        self.assertClose(quantiles[3,:], numpy.percentile(cat['parameters'], 80, axis=0), rtol=2E-3)

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(SampleSetTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
