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
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetection
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as geomEllipses
import lsst.meas.multifit as mf
import lsst.meas.multifit.utils as mfUtils
import lsst.afw.math.shapelets as shapelets

import lsst.utils.tests as utilsTests
import numpy
import unittest
import os, sys

class CompoundShapeletModelBasisTest(unittest.TestCase):
    def testIO(self):
        components = mf.CompoundShapeletBuilder.ComponentVector()
        for i in range(5):
            components.push_back(mf.ShapeletModelBasis.make(i))

        builder =  mf.CompoundShapeletBuilder(components)     
        forward = numpy.zeros_like(builder.getForward())
        reverse = numpy.zeros_like(builder.getReverse())
        i = 0
        shape=forward.shape
        area = shape[0]*shape[1]
        for r in range(shape[0]):
            for c in range(shape[1]):
                forward[r,c]=i
                reverse[r,c]=area-i
                i+=1
        builder.setMapping(forward, reverse)
        saver = builder.build()

        print forward.shape
        print reverse.shape

        filename = os.path.join("tests", "compound_shapelet.boost")
        saver.save(filename)

        loader = mf.CompoundShapeletModelBasis.load(filename)
        self.assertEqual(loader.getSize(), saver.getSize())
        loadComponents = loader.extractComponents()

        self.assertEqual(loadComponents.size(), components.size())
        for (i, j) in zip(loadComponents, components):
            self.assertAlmostEqual(i.getScale(), j.getScale())
            self.assertEqual(i.getOrder(), j.getOrder())            

        loadForward = loader.getForward()
        loadReverse = loader.getReverse()

        for r in range(shape[0]):
            for c in range(shape[1]):
                self.assertAlmostEqual(loadReverse[r,c], reverse[r,c])
                self.assertAlmostEqual(loadForward[r,c], forward[r,c])

        os.remove(filename)

    def testIntegration(self):
        components = mf.CompoundShapeletBuilder.ComponentVector()
        components.push_back(mf.ShapeletModelBasis.make(4, 1.0))
        components.push_back(mf.ShapeletModelBasis.make(4, 1.1))
        builder = mf.CompoundShapeletBuilder(components)
        builder.orthogonalize()
        basis = builder.build()
        ellipse = geomEllipses.Ellipse(geomEllipses.Axes(10.0, 8.0, 0.3), afwGeom.Point2D(0.0, 0.0))
        bounds = geomEllipses.Ellipse(ellipse)
        bounds.scale(8.0)
        footprint = afwDetection.Footprint(bounds)
        matrix = numpy.zeros((footprint.getArea(), basis.getSize()), dtype=float)
        basis.evaluate(matrix, footprint, ellipse)
        integration = numpy.zeros(basis.getSize(), dtype=float)
        basis.integrate(integration)
        factor = ellipse.getCore().getArea() / numpy.pi
        total = matrix.sum(axis=0) / factor
        self.assert_(numpy.allclose(integration, total))
        

def suite():
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(CompoundShapeletModelBasisTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
