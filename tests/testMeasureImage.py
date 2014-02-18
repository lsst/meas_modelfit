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
import matplotlib
import copy

import lsst.pex.logging
import lsst.utils.tests
import lsst.shapelet.tests
import lsst.afw.geom.ellipses
import lsst.afw.image
import lsst.afw.detection
import lsst.meas.multifit
import lsst.meas.multifit.display
import lsst.afw.display.ds9
from lsst.meas.extensions.multiShapelet import FitPsfAlgorithm

numpy.random.seed(500)

# Set to 7 for per-object messages, 10 for per-sample
lsst.pex.logging.Debug("meas.multifit.AdaptiveImportanceSampler", 0)
lsst.pex.logging.Debug("meas.multifit.TruncatedGaussian", 0)
lsst.pex.logging.Debug("meas.multifit.optimizer", 0)

DO_MAKE_PLOTS = True

DATA_DIR = os.path.join(os.environ["MEAS_MULTIFIT_DIR"], "tests", "data")

class FakeDataRef(object):

    def __init__(self, tag=None):
        self.data = dict()
        self.data['refcat'] = lsst.afw.table.SimpleCatalog.readFits(os.path.join(DATA_DIR, 'refcat.fits'))
        self.data['src'] = lsst.afw.table.SourceCatalog.readFits(os.path.join(DATA_DIR, 'src.fits'))
        self.data['calexp'] = lsst.afw.image.ExposureF(os.path.join(DATA_DIR, 'calexp.fits'))
        self.tag = tag

    def get(self, name, tag=None, immediate=True):
        tag = tag if tag is not None else self.tag
        if tag is None:
            return self.data[name]
        else:
            return self.data[name, tag]

    def put(self, obj, name, tag=None):
        tag = tag if tag is not None else self.tag
        if tag is None:
            self.data[name] = obj
        else:
            self.data[name, tag] = obj

    def datasetExists(self, name, tag=None):
        tag = tag if tag is not None else self.tag
        if tag is None:
            return name in self.data
        else:
            return (name, tag) in self.data

    def dataRef(self, name, tag=None):
        r = FakeDataRef(tag=tag)
        r.data = self.data
        return r

class MeasureImageTestCase(lsst.shapelet.tests.ShapeletTestCase):

    def setUp(self):
        self.dataRef = FakeDataRef()
        self.config = lsst.meas.multifit.MeasureCcdTask.ConfigClass()
        self.config.progressChunk = 1
        self.config.doRaise = True
        self.models = [
            'bulge+disk',
            'fixed-sersic',
            ]

    def testSampler(self):
        self.config.fitter.retarget(lsst.meas.multifit.AdaptiveImportanceSamplerTask)
        for model in self.models:
            self.config.model.name = model
            config1 = copy.deepcopy(self.config)
            config1.tag = "marginal+%s" % model
            config1.fitter.doMarginalizeAmplitudes = True
            config1.freeze()
            task1 = lsst.meas.multifit.MeasureCcdTask(config=config1, butler=self.dataRef,
                                                      name=('testMarginalSampler/%s' % model))
            task1.writeSchemas(butler=self.dataRef)
            task1.writeConfig(butler=self.dataRef)
            results1 = task1.run(self.dataRef)
            for outRecord in results1.outCat:
                self.assert_(numpy.isfinite(outRecord['fit.nonlinear']).all())
                if False:  # not yet implemented, but we should enable this test someday
                    self.assert_(numpy.isfinite(outRecord['fit.amplitudes']).all())
            if task1.model.getAmplitudeDim() > 1:
                # Direct sampling doesn't yet handle the degeneracies that can arise with
                # multi-component models very well.
                continue
            config2 = copy.deepcopy(self.config)
            config2.fitter.doMarginalizeAmplitudes = False
            config2.fitter.maxRetries = 2  # TODO: investigate why this fails with maxRetries=0
            config2.previous = config1.tag
            config2.tag = "direct+%s" % model
            config2.freeze()
            task2 = lsst.meas.multifit.MeasureCcdTask(config=config2, butler=self.dataRef,
                                                      name=('testDirectSampler/%s' % model))
            task2.writeSchemas(butler=self.dataRef)
            task2.writeConfig(butler=self.dataRef)
            results2 = task2.run(self.dataRef)
            for outRecord in results2.outCat:
                self.assert_(numpy.isfinite(outRecord['fit.nonlinear']).all())
                self.assert_(numpy.isfinite(outRecord['fit.amplitudes']).all())

    def testOptimizer(self):
        self.config.fitter.retarget(lsst.meas.multifit.OptimizerTask)
        self.config.fitter.doRecordHistory = True
        for model in self.models:
            self.config.model.name = model
            name = 'testOptimizer/%s' % model
            task = lsst.meas.multifit.MeasureCcdTask(config=self.config, name=name)
            results = task.run(self.dataRef)
            for outRecord in results.outCat:
                self.assertFalse(outRecord.get("fit.flags"))
                self.assert_(numpy.isfinite(outRecord['fit.nonlinear']).all())
                self.assert_(numpy.isfinite(outRecord['fit.amplitudes']).all())

    def tearDown(self):
        del self.dataRef
        del self.config

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
