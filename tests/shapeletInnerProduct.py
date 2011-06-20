#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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
import lsst.afw.math.shapelets
import lsst.afw.detection
import lsst.afw.geom.ellipses
import lsst.meas.multifit.utils
import lsst.utils.tests
import numpy
import unittest

# NOTE: this is a regression test that depends on some floating point-comparisons
#       that don't always succeed with different random numbers.  If it fails,
#       look at the example with the same name to verify if this is really a
#       failure or just an overly optimistic comparison between and analytical
#       and numerical computations.
numpy.random.seed(1)

class ShapeletInnerProductTestCase(unittest.TestCase):

    def assertClose(self, a, b, rtol=1E-5, atol=1E-8):
        self.assert_(numpy.allclose(a, b, rtol=rtol, atol=atol), "%s\n != \n%s" % (a, b))

    def compare(self, basis, ellipse, factor=15):
        n = basis.getSize()
        bounds = lsst.afw.geom.ellipses.Ellipse(ellipse)
        bounds.scale(factor)
        box = lsst.afw.geom.Box2I(bounds.computeEnvelope())
        footprint = lsst.afw.detection.Footprint(box)
        array = numpy.zeros((footprint.getArea(), basis.getSize()), dtype=float)
        basis.evaluate(array, footprint, ellipse)
        shape = (box.getHeight(), box.getWidth())
        m1 = basis.computeInnerProductMatrix(ellipse.getCore())
        m2 = numpy.zeros_like(m1)
        for i in range(n):
            for j in range(i + 1):
                image = (array[:,i] * array[:,j]).reshape(*shape)
                m2[i,j] = image.sum()
                m2[j,i] = m2[i,j]
        mask = numpy.abs(m1) > 1
        d = numpy.abs((m1[mask] - m2[mask]) / (m1[mask]))
        self.assert_(d.max() < 0.1, "%f >= %f" % (d.max(), 0.1))
    
    def testInnerProducts(self):
        basis1 = lsst.meas.multifit.ShapeletModelBasis.make(2, 1.2)
        basis2 = lsst.meas.multifit.utils.loadBasis("ed+06:2000")
        psf = lsst.afw.math.shapelets.ShapeletFunction(
            2, lsst.afw.math.shapelets.HERMITE,
            lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes(1.0, 1.2, 0.0))
            )
        psf.getCoefficients()[:] = numpy.random.randn(psf.getCoefficients().size)
        convolved1 = basis1.convolve(psf)
        convolved2 = basis2.convolve(psf)
        ellipse = lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes(7.0, 5.0, 0.0))
        self.compare(basis1, ellipse)
        self.compare(basis2, ellipse)
        self.compare(convolved1, ellipse)
        self.compare(convolved2, ellipse)


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """Returns a suite containing all the test cases in this module."""

    lsst.utils.tests.init()

    suites = []
    suites += unittest.makeSuite(ShapeletInnerProductTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    lsst.utils.tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
