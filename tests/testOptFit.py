#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2015 LSST/AURA
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

import lsst.afw.geom
import lsst.afw.table
import lsst.utils.tests
import lsst.meas.modelfit

from lsst.meas.base.tests import AlgorithmTestCase, TestDataset

class OptFitTestCase(AlgorithmTestCase):
    """Test case for the OptFit measurement plugin
    """

    def setUp(self):
        AlgorithmTestCase.setUp(self)
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(0, 0),
                                        lsst.afw.geom.Extent2I(200, 100))
        self.dataset = TestDataset(self.bbox)
        # first source is a point
        self.dataset.addSource(100000.0, lsst.afw.geom.Point2D(50.1, 49.8))
        # second source is extended
        self.dataset.addSource(100000.0, lsst.afw.geom.Point2D(149.9, 50.3),
                               lsst.afw.geom.ellipses.Quadrupole(8, 9, 3))

    def tearDown(self):
        AlgorithmTestCase.tearDown(self)
        del self.bbox
        del self.dataset

    def testPlugins(self):
        plugin = "modelfit_OptFit"
        dependencies = ("modelfit_ShapeletPsfApprox",)
        sfmConfig = self.makeSingleFrameMeasurementConfig(plugin, dependencies=dependencies)
        sfmConfig.plugins["modelfit_ShapeletPsfApprox"].sequence.append("Full")
        sfmTask = self.makeSingleFrameMeasurementTask(config=sfmConfig)
        # catalog1 will contain both the SFM outputs and the truth catalog for sources in exposure 1.
        # Those SFM outputs will also be used as the references for the forced task.
        exposure, catalog = self.dataset.realize(10.0, sfmTask.schema)
        sfmTask.run(catalog, exposure)

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(OptFitTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
