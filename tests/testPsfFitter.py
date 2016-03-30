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
import glob
import numpy

import lsst.utils.tests
import lsst.shapelet
import lsst.afw.geom.ellipses
import lsst.meas.modelfit
import lsst.meas.base

numpy.random.seed(500)

lsst.pex.logging.Debug("meas.modelfit.PsfFitter", 10)

DATA_DIR = os.path.join(os.environ["MEAS_MODELFIT_DIR"], "tests", "data")

def computeMoments(image):
    """Helper function to compute moments of a postage stamp about its origin."""
    maskedImage = lsst.afw.image.MaskedImageD(image)
    result = lsst.meas.base.SdssShapeAlgorithm.computeAdaptiveMoments(
        maskedImage,
        lsst.afw.geom.Point2D(0.0, 0.0)
        )
    return result.getShape()

class PsfFitterTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.configs = {}
        self.configs['ellipse'] = lsst.meas.modelfit.PsfFitterConfig()
        self.configs['full'] = lsst.meas.modelfit.PsfFitterConfig()
        self.configs['full'].inner.order = 0
        self.configs['full'].primary.order = 4
        self.configs['full'].wings.order = 4
        self.configs['full'].outer.order = 0

    def tearDown(self):
        del self.configs

    def testApply(self):
        tolerances = {"full": 3E-4, "ellipse": 8E-3}
        from matplotlib import pyplot
        import matplotlib.colors
        for filename in glob.glob(os.path.join(DATA_DIR, "psfs", "great3*.fits")):
            kernelImage = lsst.afw.image.ImageD(filename)
            shape = computeMoments(kernelImage)
            vrange = kernelImage.getArray().max() - kernelImage.getArray().min()
            vmax = kernelImage.getArray().max() + vrange*0.05
            vmin = kernelImage.getArray().min() - vrange*0.05
            norm = matplotlib.colors.SymLogNorm(vrange*0.1, vmin=-vmax, vmax=vmax)
            for configKey in ["full", "ellipse"]:
                fitter = lsst.meas.modelfit.PsfFitter(self.configs[configKey].makeControl())
                multiShapeletFit = fitter.makeInitial(shape)
                fitter.apply(multiShapeletFit, kernelImage)
                pyplot.figure("%s: %s" % (filename, configKey))
                dump = lsst.meas.modelfit.debugRegistry["psf.fitEM"]
                nIterations = len(dump)/3
                nComponents = fitter.getComponentCount()
                for iteration in xrange(nIterations):
                    matrix = dump["iter=%05d:matrix" % iteration].transpose().reshape(
                        nComponents,
                        kernelImage.getHeight(),
                        kernelImage.getWidth()
                    )
                    solution = dump["iter=%05d:solution" % iteration]
                    model = dump["iter=%05d:model" % iteration].reshape(
                        kernelImage.getHeight(),
                        kernelImage.getWidth()
                    )
                    for n in xrange(nComponents):
                        ax = pyplot.subplot(nComponents+2, nIterations, n*nIterations + iteration + 1)
                        ax.imshow(matrix[n]*solution[n], origin='lower', interpolation='nearest',
                                  norm=norm)
                        ax.axis("off")
                        if n == 0:
                            ax.set_title(iteration)
                    ax = pyplot.subplot(nComponents+2, nIterations,
                                        nComponents*nIterations + iteration + 1)
                    ax.imshow(model, origin='lower', interpolation='nearest', norm=norm)
                    ax.axis("off")
                    ax = pyplot.subplot(nComponents+2, nIterations,
                                        (nComponents + 1)*nIterations + iteration + 1)
                    ax.imshow(model - kernelImage.getArray(), origin='lower', interpolation='nearest',
                              norm=norm)
                    ax.axis("off")
                pyplot.show()
                modelImage = lsst.afw.image.ImageD(kernelImage.getBBox(lsst.afw.image.PARENT))
                multiShapeletFit.evaluate().addToImage(modelImage)
                self.assertClose(kernelImage.getArray(), modelImage.getArray(),
                                 atol=tolerances[configKey],
                                 plotOnFailure=False)

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(PsfFitterTestCase)
    suites += unittest.makeSuite(lsst.utils.tests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
