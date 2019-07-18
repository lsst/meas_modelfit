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
import os
import numpy

import lsst.utils.tests
import lsst.meas.modelfit


class IntegralsTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        numpy.random.seed(500)

    def testBVN(self):
        data = numpy.loadtxt(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                          "reference", "bvn.txt"), delimiter=',')
        for h, k, r, p1 in data:
            p2 = lsst.meas.modelfit.detail.bvnu(h, k, r)
            self.assertFloatsAlmostEqual(p1, p2, rtol=1E-14)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
