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

class ObjectiveTestCase(lsst.shapelet.tests.ShapeletTestCase):

    def testZeroRadius(self):
        psf = self.makeRandomMultiShapeletFunction()
        psf.normalize()
        bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-50,50), lsst.afw.geom.Point2I(50,50))
        image = lsst.afw.image.MaskedImageF(bbox)
        image.getImage().getArray()[:,:] = numpy.random.randn(bbox.getHeight(), bbox.getWidth())
        image.getVariance().getArray()[:,:] = 1.0
        footprint = lsst.afw.detection.Footprint(bbox)
        ModelConfig = lsst.meas.multifit.models.BulgeDiskModelConfig
        basis = ModelConfig.makeBasis(ModelConfig())
        objective = lsst.meas.multifit.SingleEpochObjective(
            lsst.meas.multifit.SingleEpochObjectiveControl(),
            basis, psf, image, footprint
            )
        ellipse = lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Quadrupole(0.0, 0.0, 0.0))
        lg = objective.evaluate(ellipse)
        self.assertFalse(numpy.isnan(lg.grad).any())
        self.assertFalse(numpy.isnan(lg.fisher).any())

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(ObjectiveTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
