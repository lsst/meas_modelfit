from builtins import zip
#
# LSST Data Management System
#
# Copyright 2008-2016  AURA/LSST.
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
# see <https://www.lsstcorp.org/LegalNotices/>.
#
import unittest

import lsst.afw.geom
import lsst.afw.table
import lsst.utils.tests
import lsst.meas.modelfit

import lsst.meas.base.tests


# n.b. Some tests here depend on the noise realization in the test data
# or from the numpy random number generator.
# For the current test data and seed value, they pass, but they may not
# if the test data is regenerated or the seed value changes.  I've marked
# these with an "rng dependent" comment.  In most cases, they test that
# the measured flux lies within 2 sigma of the correct value, which we
# should expect to fail sometimes.

class CModelTestCase(lsst.meas.base.tests.AlgorithmTestCase, lsst.utils.tests.TestCase):
    """Test case for the CModel measurement plugins
    """

    def setUp(self):
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(0, 0),
                                        lsst.afw.geom.Extent2I(200, 100), invert=False)
        self.dataset = lsst.meas.base.tests.TestDataset(self.bbox)
        # first source is a point
        self.dataset.addSource(100000.0, lsst.afw.geom.Point2D(50.1, 49.8))
        # second source is extended
        self.dataset.addSource(100000.0, lsst.afw.geom.Point2D(149.9, 50.3),
                               lsst.afw.geom.ellipses.Quadrupole(8, 9, 3))

    def tearDown(self):
        del self.bbox
        del self.dataset

    def checkOutputs(self, measCat, truthCat=None):
        """Test that the outputs of the CModel plugins are reasonable, and that the bookkeeping works.

        Science-quality tests either in testCModel.py (where we call the same code via a different interface)
        or something we have to do statistically on real data.
        """
        if truthCat is None:
            truthCat = measCat
        for measRecord, truthRecord in zip(measCat, truthCat):
            trueFlux = truthRecord.get("truth_flux")
            self.assertFalse(measRecord.get("modelfit_CModel_initial_flag"))
            self.assertFalse(measRecord.get("modelfit_CModel_exp_flag"))
            self.assertFalse(measRecord.get("modelfit_CModel_dev_flag"))
            self.assertFalse(measRecord.get("modelfit_CModel_flag"))
            self.assertFloatsAlmostEqual(measRecord.get("modelfit_CModel_flux"), trueFlux, rtol=0.5)
            self.assertGreater(measRecord.get("modelfit_CModel_fluxSigma"), 0.0)
            self.assertFloatsAlmostEqual(measRecord.get("modelfit_CModel_initial_flux"), trueFlux, rtol=0.5)
            self.assertGreater(measRecord.get("modelfit_CModel_initial_fluxSigma"), 0.0)
            self.assertFloatsAlmostEqual(measRecord.get("modelfit_CModel_exp_flux"), trueFlux, rtol=0.5)
            self.assertGreater(measRecord.get("modelfit_CModel_exp_fluxSigma"), 0.0)
            self.assertFloatsAlmostEqual(measRecord.get("modelfit_CModel_dev_flux"), trueFlux, rtol=0.5)
            self.assertGreater(measRecord.get("modelfit_CModel_dev_fluxSigma"), 0.0)

    def testPlugins(self):
        """Test that the plugin for single-frame measurement works, then use those outputs to test
        that the forced measurement plugin works."""
        plugin = "modelfit_CModel"
        dependencies = ("modelfit_DoubleShapeletPsfApprox", "base_PsfFlux")
        sfmTask = self.makeSingleFrameMeasurementTask(plugin, dependencies=dependencies)
        forcedTask = self.makeForcedMeasurementTask(plugin, dependencies=dependencies,
                                                    refSchema=sfmTask.schema)
        # catalog1 will contain both the SFM outputs and the truth catalog for sources in exposure 1.
        # Those SFM outputs will also be used as the references for the forced task.
        exposure1, catalog1 = self.dataset.realize(10.0, sfmTask.schema, randomSeed=0)
        sfmTask.run(catalog1, exposure1)
        self.checkOutputs(catalog1)
        if False:  # this line should be re-enabled on DM-5405
            wcs2 = self.dataset.makePerturbedWcs(self.dataset.exposure.getWcs(), randomSeed=0)
        else:
            wcs2 = self.dataset.exposure.getWcs()
        dataset2 = self.dataset.transform(wcs2)
        # catalog2 will contain only the truth catalog for sources in exposure 1; the structure of
        # ForcedMeasurementTask means we can't put the outputs in the same catalog.
        exposure2, catalog2 = dataset2.realize(10.0, dataset2.makeMinimalSchema(), randomSeed=1)
        refWcs = exposure1.getWcs()
        refCat = catalog1
        measCat = forcedTask.generateMeasCat(exposure2, refCat, refWcs)
        forcedTask.attachTransformedFootprints(measCat, refCat, exposure2, refWcs)
        forcedTask.run(measCat, exposure2, refCat, refWcs)
        self.checkOutputs(measCat, catalog2)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
