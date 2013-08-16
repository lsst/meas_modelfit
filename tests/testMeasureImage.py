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

import os
import unittest
import numpy

import lsst.utils.tests
import lsst.afw.geom.ellipses
import lsst.afw.image
import lsst.afw.detection
import lsst.shapelet.tests
import lsst.meas.multifit

numpy.random.seed(500)

DATA_DIR = os.path.join(os.environ["MEAS_MULTIFIT_DIR"], "tests", "data")

class MeasureImageTestCase(lsst.shapelet.tests.ShapeletTestCase):

    def setUp(self):
        self.task = lsst.meas.multifit.MeasureCcdTask()
        self.modelfits = lsst.meas.multifit.ModelFitCatalog.readFits(os.path.join(DATA_DIR, "outcat.fits"))
        self.calexp = lsst.afw.image.ExposureF(os.path.join(DATA_DIR, "calexp.fits"))

    def tearDown(self):
        del self.task
        del self.modelfits
        del self.calexp

    def testProcessObject(self):
        # For now, just test that there are no exceptions; should test the outputs better later
        result = self.task.processObject(exposure=self.calexp, record=self.modelfits[0])

        # Test persistence of ModelFitCatalog and especially SampleSet (done here just because
        # it's otherwise a pain to build a realistically complex SampleSet).
        filename = "testModelFitPersistence.fits"
        self.modelfits.writeFits(filename)
        loaded = lsst.meas.multifit.ModelFitCatalog.readFits(filename)
        self.assertEqual(self.modelfits.schema.compare(loaded.schema, lsst.afw.table.Schema.IDENTICAL),
                         lsst.afw.table.Schema.IDENTICAL)
        self.assertEqual(len(self.modelfits), len(loaded))
        samples1 = self.modelfits[0].getSamples()
        samples2 = loaded[0].getSamples()
        cat1 = samples1.getCatalog().copy(deep=True)
        cat2 = samples2.getCatalog().copy(deep=True)
        self.assertEqual(len(cat1), len(cat2))
        self.assertEqual(samples1.getParameterDefinition(), samples2.getParameterDefinition())
        self.assertEqual(samples1.getDataSquaredNorm(), samples2.getDataSquaredNorm())
        # n.b. just using assertClose because it lets us test arrays
        self.assertClose(cat1.get("joint.grad"), cat2.get("joint.grad"), rtol=0.0, atol=0.0)
        self.assertClose(cat1.get("marginal"), cat2.get("marginal"), rtol=0.0, atol=0.0)
        self.assertClose(cat1.get("proposal"), cat2.get("proposal"), rtol=0.0, atol=0.0)
        self.assertClose(cat1.get("weight"), cat2.get("weight"), rtol=0.0, atol=0.0)
        self.assertClose(cat1.get("parameters"), cat2.get("parameters"), rtol=0.0, atol=0.0)
        self.assertClose(cat1.get("joint.fisher"), cat2.get("joint.fisher"), rtol=0.0, atol=0.0)
        os.remove(filename)

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(MeasureImageTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
