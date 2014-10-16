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

import lsst.afw.geom
import lsst.afw.table
import lsst.utils.tests
import lsst.meas.base.tests
import lsst.meas.modelfit

import lsst.afw.display

numpy.random.seed(1234567)

# n.b. Some tests here depend on the noise realization in the test data
# or from the numpy random number generator.
# For the current test data and seed value, they pass, but they may not
# if the test data is regenerated or the seed value changes.  I've marked
# these with an "rng dependent" comment.  In most cases, they test that
# the measured flux lies within 2 sigma of the correct value, which we
# should expect to fail sometimes.

class CModelTestCase(lsst.meas.base.tests.AlgorithmTestCase):
    """Test case for the CModel measurement plugins
    """

    def setUp(self):
        lsst.meas.base.tests.AlgorithmTestCase.setUp(self)
        self.record = self.truth[0]

    def tearDown(self):
        lsst.meas.base.tests.AlgorithmTestCase.tearDown(self)
        del self.record

    def checkOutputs(self, measCat):
        """Test that the outputs of the CModel plugins are reasonable, and that the bookkeeping works.

        Science-quality tests either in testCModel.py (where we call the same code via a different interface)
        or something we have to do statistically on real data.
        """
        for measRecord, truthRecord in zip(measCat, self.truth):
            trueFlux = truthRecord.get("truth_flux")
            if not measRecord.getShapeFlag():
                self.assertFalse(measRecord.get("modelfit_CModel_initial_flag"))
                self.assertFalse(measRecord.get("modelfit_CModel_exp_flag"))
                self.assertFalse(measRecord.get("modelfit_CModel_dev_flag"))
                self.assertFalse(measRecord.get("modelfit_CModel_flag"))
                self.assertClose(measRecord.get("modelfit_CModel_flux"), trueFlux, rtol=0.5)
                self.assertGreater(measRecord.get("modelfit_CModel_fluxSigma"), 0.0)
                self.assertClose(measRecord.get("modelfit_CModel_initial_flux"), trueFlux, rtol=0.5)
                self.assertGreater(measRecord.get("modelfit_CModel_initial_fluxSigma"), 0.0)
                self.assertClose(measRecord.get("modelfit_CModel_exp_flux"), trueFlux, rtol=0.5)
                self.assertGreater(measRecord.get("modelfit_CModel_exp_fluxSigma"), 0.0)
                self.assertClose(measRecord.get("modelfit_CModel_dev_flux"), trueFlux, rtol=0.5)
                self.assertGreater(measRecord.get("modelfit_CModel_dev_fluxSigma"), 0.0)
            else:
                self.assertTrue(measRecord.get("modelfit_CModel_initial_flag"))
                self.assertTrue(measRecord.get("modelfit_CModel_exp_flag"))
                self.assertTrue(measRecord.get("modelfit_CModel_dev_flag"))
                self.assertTrue(measRecord.get("modelfit_CModel_flag"))

    def testPlugins(self):
        # Start with a run on some simulated data, using the single-frame measurement driver
        sfmConfig = lsst.meas.base.SingleFrameMeasurementTask.ConfigClass()
        sfmConfig.plugins.names = ["base_SdssShape", "base_PsfFlux", "modelfit_CModel",
                                   "modelfit_ShapeletPsfApprox"]
        sfmConfig.slots.centroid = None
        sfmConfig.slots.shape = "base_SdssShape"
        sfmConfig.slots.psfFlux = "base_PsfFlux"
        sfmConfig.slots.apFlux = None
        sfmConfig.slots.modelFlux = None
        sfmConfig.slots.instFlux = None
        sfmSchemaMapper = lsst.afw.table.SchemaMapper(self.truth.schema)
        sfmSchemaMapper.addMinimalSchema(self.truth.schema)
        sfmTask = lsst.meas.base.SingleFrameMeasurementTask(config=sfmConfig,
                                                            schema=sfmSchemaMapper.editOutputSchema())
        sfmMeasCat = lsst.afw.table.SourceCatalog(sfmSchemaMapper.getOutputSchema())
        sfmMeasCat.extend(self.truth, sfmSchemaMapper)
        sfmMeasCat.table.defineCentroid("truth")
        sfmTask.run(sfmMeasCat, self.calexp)
        self.checkOutputs(sfmMeasCat)

        # Now we use the SFM results as the reference catalog for a forced measurement run
        forcedConfig = lsst.meas.base.ForcedMeasurementTask.ConfigClass()
        forcedConfig.plugins.names = ["base_TransformedCentroid", "base_TransformedShape",
                                      "base_PsfFlux", "modelfit_CModel", "modelfit_ShapeletPsfApprox"]
        forcedConfig.slots.centroid = "base_TransformedCentroid"
        forcedConfig.slots.shape = "base_TransformedShape"
        forcedConfig.slots.psfFlux = "base_PsfFlux"
        forcedConfig.slots.apFlux = None
        forcedConfig.slots.modelFlux = None
        forcedConfig.slots.instFlux = None
        forcedTask = lsst.meas.base.ForcedMeasurementTask(config=forcedConfig, refSchema=sfmMeasCat.schema)
        forcedMeasCat = forcedTask.run(self.calexp, sfmMeasCat, self.calexp.getWcs()).sources
        # only test the first three simulated sources: the last two are deblend children, and because
        # of the lack of heavy footprints in forced mode, we don't measure those well at all.
        self.checkOutputs(forcedMeasCat[:3])

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(CModelTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
