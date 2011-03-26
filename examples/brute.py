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
import lsst.afw.image
import lsst.afw.geom.ellipses
import lsst.meas.multifit.viewer
import lsst.meas.multifit.sampling
import numpy
from matplotlib import pyplot

def main(xr, yr, sn=20.0):
    viewer = lsst.meas.multifit.viewer.StarViewer.makeExample(sn=sn)
    grid = numpy.zeros(yr.shape + xr.shape, dtype=float)
    for iy, y in enumerate(yr):
        for ix, x in enumerate(xr):
            parameters = numpy.array([x, y], dtype=float)
            grid[iy, ix] = viewer.update(parameters)
    return grid
            
def student(self, x, dof, mu, sigma_inv):
    y = numpy.dot(sigma_inv, x - mu)
    z = numpy.sum(y**2, axis=1)
    return (1 + z / dof)**(-0.5*(dof + mu.size))
