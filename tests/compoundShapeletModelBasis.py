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
        mapping = numpy.zeros_like(builder.getMapping())
        i = 0
        shape=mapping.shape
        for r in range(shape[0]):
            for c in range(shape[1]):
                mapping[r,c]=i
                i+=1
        builder.setMapping(mapping)
        saver = builder.build()

        print mapping.shape

        filename = os.path.join("tests", "compound_shapelet.boost")
        saver.save(filename)

        loader = mf.CompoundShapeletModelBasis.load(filename)
        self.assertEqual(loader.getSize(), saver.getSize())
        loadComponents = loader.extractComponents()

        self.assertEqual(loadComponents.size(), components.size())
        for (i, j) in zip(loadComponents, components):
            self.assertAlmostEqual(i.getScale(), j.getScale())
            self.assertEqual(i.getOrder(), j.getOrder())            

        loadMapping = loader.getMapping()

        for r in range(shape[0]):
            for c in range(shape[1]):
                self.assertAlmostEqual(loadMapping[r,c], mapping[r,c])

        os.remove(filename)

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
