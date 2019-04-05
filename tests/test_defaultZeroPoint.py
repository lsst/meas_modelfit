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
import numpy as np

import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.meas.modelfit as measModel
import lsst.utils


class DefaultZeroPointTestCase(lsst.utils.tests.TestCase):
    '''Test that photoCalib objects containing invalid calibrations are replaced with a default having a
    zero-point of magnitude 27.
    '''

    def setUp(self):
        '''Create two calibs, one with a valid zero-point and one without. Use these to create two UnitSystem
        objects.
        '''
        self.mag2Flux = lambda m: 10.0**(m/2.5)
        self.flux2Mag = lambda f: 2.5*np.log10(f)

        photoCalibNoZero = afwImage.PhotoCalib()
        photoCalibWithZero = afwImage.makePhotoCalibFromCalibZeroPoint(self.mag2Flux(25))

        scale = 0.2 * afwGeom.arcseconds
        wcs = afwGeom.makeSkyWcs(crpix=afwGeom.Point2D(),
                                 crval=afwGeom.SpherePoint(45.0, 45.0, afwGeom.degrees),
                                 cdMatrix=afwGeom.makeCdMatrix(scale=scale))

        self.unitNoZero = measModel.UnitSystem(wcs, photoCalibNoZero)
        self.unitWithZero = measModel.UnitSystem(wcs, photoCalibWithZero)

    def testZeroPoint(self):
        '''Verify that the UnitSystem with an invalid photoCalib gets populated with a valid photoCalib, and that
        nothing happens to the UnitSystem with a valid photoCalib.
        '''
        magNoZero = self.flux2Mag(self.unitNoZero.photoCalib.getInstFluxAtZeroMagnitude())
        magWithZero = self.flux2Mag(self.unitWithZero.photoCalib.getInstFluxAtZeroMagnitude())

        print((magNoZero, magWithZero))

        self.assertAlmostEqual(magNoZero, 27, 6)
        self.assertAlmostEqual(magWithZero, 25, 6)

    def tearDown(self):
        del self.unitNoZero
        del self.unitWithZero


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
