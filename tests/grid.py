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
import lsst.afw.image
import lsst.afw.detection
import lsst.afw.geom
import lsst.afw.geom.ellipses
import lsst.meas.multifit as mf

import lsst.utils.tests as utilsTests
import numpy
import unittest

class DefinitionTest(unittest.TestCase):

    def setUp(self):
        self.galaxyEllipse = lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes(4.0, 3.0, 0.25), 
                                                            lsst.afw.geom.Point2D(10.3, 42.1))
        self.galaxy = mf.definition.Object.makeGalaxy(
            1,
            mf.ShapeletModelBasis.make(2),
            self.galaxyEllipse,
            True, True, True
            )
        self.starPoint = lsst.afw.geom.Point2D(12.1, 56.7)
        self.star = mf.definition.Object.makeStar(3, self.starPoint, False, True)
        self.footprint = lsst.afw.detection.Footprint(self.galaxyEllipse)
        self.data = numpy.random.randn(self.footprint.getArea())
        self.frame = mf.definition.Frame(0, self.footprint, self.data)

    def tearDown(self):
        del self.frame
        del self.footprint

    def testComponents(self):
        position1 = mf.Position(4.3, 2.0)
        position2 = mf.Position(-2.1, 0.5)
        p = mf.definition.PositionComponent.make(position1)
        self.assert_(p.isActive())
        self.assertEqual(p.getValue(), position1)
        p.setValue(position2)
        self.assertEqual(p.getValue(), position2)
        p.getValue().setX(3.0)
        self.assertEqual(p.getValue().getX(), 3.0)
        p.setActive(False)
        self.assertFalse(p.isActive())
        p.getBounds().max = 3.0
        self.assertEqual(p.getBounds().max, 3.0)

        ellipticity1 = mf.Ellipticity(0.6, -0.1)
        ellipticity2 = mf.Ellipticity(-0.25, 0.3)
        e = mf.definition.EllipticityComponent.make(ellipticity1, False)
        self.assertEqual(e.getValue().getE1(), ellipticity1.getE1())
        self.assertEqual(e.getValue().getE2(), ellipticity1.getE2())
        self.assertFalse(e.isActive())
        e.setValue(ellipticity2)
        self.assertEqual(e.getValue().getE1(), ellipticity2.getE1())
        self.assertEqual(e.getValue().getE2(), ellipticity2.getE2())
        e.getValue().setE1(0.4)
        self.assertEqual(e.getValue().getE1(), 0.4)
        e.setActive(True)
        self.assert_(e.isActive())
        
        radius1 = mf.Radius(2.3)
        radius2 = mf.Radius(1.6)
        r = mf.definition.RadiusComponent.make(radius1, True)
        self.assertEqual(float(r.getValue()), float(radius1))
        self.assert_(r.isActive())
        r.setValue(radius2)
        self.assertEqual(float(r.getValue()), float(radius2))
        r.setActive(False)
        self.assertFalse(r.isActive())
        r.setBounds(r.Bounds(0.0, 5.0))
        self.assertEqual(r.getBounds().min, 0.0)
        self.assertEqual(r.getBounds().max, 5.0)

    def testParameterSharing(self):
        self.assertEqual(self.star.getPosition().getValue(), self.starPoint)
        self.assertEqual(self.galaxy.getPosition().getValue(), self.galaxyEllipse.getCenter())
        self.galaxy.setPosition(self.star.getPosition())
        self.assertEqual(self.galaxy.getPosition(), self.star.getPosition())
        self.assertNotEqual(self.galaxy.getPosition(), self.star.getRadius()) # this shouldn't throw

    def testObjectAccessors(self):
        self.assertEqual(self.star.id, 3)
        self.assertEqual(self.galaxy.id, 1)
        self.assertFalse(self.star.isVariable())
        self.star.setVariable(True)
        self.assert_(self.star.isVariable())
        self.assertEqual(self.star.getBasis(), None)
        self.star.setBasis(mf.ShapeletModelBasis.make(0))
        self.assertEqual(self.star.getBasis().getSize(), 1)
        self.assertEqual(self.galaxy.getRadiusFactor(), 1.0)
        self.galaxy.setRadiusFactor(2.0)
        self.assertEqual(self.galaxy.getRadiusFactor(), 2.0)

    def testObjectSet(self):
        s = mf.definition.ObjectSet()
        self.assert_(s.insert(self.star))
        self.assert_(s.insert(self.galaxy))
        self.assertFalse(s.insert(self.star))
        self.assertEqual(len(s), 2)
        self.assertEqual(list(s), s.values())

    def testFrame(self):
        self.assertEqual(self.frame.getFootprint().getArea(), self.footprint.getArea())
        self.assertEqual(self.frame.getFootprint().getBBox(), self.footprint.getBBox())
        self.assert_(self.frame.getWcs() is None)
        self.assert_(self.frame.getPsf() is None)
        self.assertEqual(self.frame.getWeights().shape, (0,))
        self.assert_((self.frame.getData() == self.data).all())
        newData = numpy.random.randn(self.footprint.getArea())
        self.frame.setData(newData)
        self.assert_((self.frame.getData() == newData).all())
        weights = numpy.random.randn(self.footprint.getArea())
        self.frame.setWeights(weights)
        self.assert_((self.frame.getWeights() == weights).all())
        psf = lsst.afw.detection.createPsf("DoubleGaussian", 19, 19, 2.0, 1.0)
        self.frame.setPsf(psf)
        self.assert_(self.frame.getPsf() is not None)

    def testFrameBuilding(self):
        bbox = lsst.afw.geom.BoxI(lsst.afw.geom.PointI(25, 40), lsst.afw.geom.ExtentI(10, 30))
        psf = lsst.afw.detection.createPsf("DoubleGaussian", 19, 19, 2.0, 1.0)
        fp = lsst.afw.detection.Footprint(bbox)
        expF = lsst.afw.image.ExposureF(bbox)
        expF.setPsf(psf)
        expD = lsst.afw.image.ExposureF(bbox)
        expD.setPsf(psf)
        frameF = mf.definition.Frame.make(4, expF, fp)
        frameD = mf.definition.Frame.make(5, expD, fp)
        self.assertEqual(frameF.id, 4)
        self.assertEqual(frameD.id, 5)

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
