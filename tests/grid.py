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

import lsst.utils.tests as utilsTests
import numpy
import unittest

class DefinitionTest(unittest.TestCase):

    def setUp(self):
        self.galaxyEllipse = geomEllipses.Ellipse(geomEllipses.Axes(4.0, 3.0, 0.25), 
                                                  afwGeom.Point2D(10.3, 42.1))
        self.galaxy = mf.definition.Object.makeGalaxy(
            1,
            mf.ShapeletModelBasis.make(2),
            self.galaxyEllipse,
            True, True, True
            )
        self.starPoint = afwGeom.Point2D(12.1, 56.7)
        self.star = mf.definition.Object.makeStar(3, self.starPoint, False, True)
        self.footprint = afwDetection.Footprint(self.galaxyEllipse);
        self.data = numpy.random.randn(self.footprint.getArea())
        self.frame = mf.definition.Frame(0, self.footprint, self.data)

    def tearDown(self):
        del self.frame
        del self.footprint
        
    def testParameterSharing(self):
        self.assertEqual(self.star.position.getReference(), self.starPoint)
        self.assertEqual(self.galaxy.position.getReference(), self.galaxyEllipse.getCenter())
        self.galaxy.position = self.star.position
        self.assertEqual(self.galaxy.position, self.star.position)
        self.assertNotEqual(self.galaxy.position, self.star.radius) # this shouldn't throw

    def testObjectSet(self):
        s = mf.definition.ObjectSet()
        self.assert_(s.insert(self.star))
        self.assert_(s.insert(self.galaxy))
        self.assertFalse(s.insert(self.star))
        self.assertEqual(len(s), 2)
        self.assertEqual(list(s), s.values())

    def testFrame(self):
        self.assertEqual(self.frame.footprint.getArea(), self.footprint.getArea())
        self.assertEqual(self.frame.footprint.getBBox(), self.footprint.getBBox())
        self.assert_(self.frame.wcs is None)
        self.assert_(self.frame.psf is None)
        self.assertEqual(self.frame.weights.shape, (0,))
        self.assert_((self.frame.data == self.data).all())
        newData = numpy.random.randn(self.footprint.getArea())
        self.frame.data = newData
        self.assert_((self.frame.data == newData).all())
        weights = numpy.random.randn(self.footprint.getArea())
        self.frame.weights = weights
        self.assert_((self.frame.weights == weights).all())
        psf = afwDetection.createPsf("DoubleGaussian", 19, 19, 2.0, 1.0)
        self.frame.psf = psf
        self.assert_(self.frame.psf is not None)

    def testDefinition(self):
        definition = mf.Definition()
        definition.objects.insert(self.star)
        definition.objects.insert(self.galaxy)
        self.assertEqual(list(definition.objects), [self.galaxy, self.star]) # sorted by id
        definition.frames.insert(self.frame)
        self.assertEqual(list(definition.frames), [self.frame])

def suite():
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(DefinitionTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
