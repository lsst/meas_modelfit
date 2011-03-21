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


class EvaluatorTest(unittest.TestCase):
    def testStarConstruction(self):        
        bbox = afwGeom.BoxI(afwGeom.PointI(25, 40), afwGeom.ExtentI(10, 30))        
        exp = afwImage.ExposureF(bbox)
        psf = afwDetection.createPsf("DoubleGaussian", 19, 19, 2.0, 1.0)
        exp.setPsf(psf)
        fp = afwDetection.Footprint(bbox)
        point = afwGeom.PointD(30.5, 55.8)
        evaluator = mf.Evaluator.make(exp, fp, point)

    def testGalaxyConstruction(self):        
        axes = geomEllipses.Axes(30, 10, 0)
        point = afwGeom.PointD(30.5, 55.8)
        ellipse = geomEllipses.Ellipse(axes, point)

        bbox = afwGeom.BoxI(ellipse.computeEnvelope())
        exp = afwImage.ExposureF(bbox)
        psf = afwDetection.createPsf("DoubleGaussian", 19, 19, 2.0, 1.0)
        exp.setPsf(psf)
        fp = afwDetection.Footprint(ellipse);
        basis = mf.ShapeletModelBasis.make(5, 1.0)
        evaluator = mf.Evaluator.make(exp, fp, basis, ellipse, False, True, True)

def suite():
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(EvaluatorTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
