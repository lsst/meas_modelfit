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

def makeGrid(viewer, ranges):
    assert(len(ranges) == viewer.evaluator.getParameterSize())
    nParameters = len(ranges)
    grid = numpy.zeros(tuple(len(r) for r in ranges), dtype=float)
    parameters = numpy.zeros(nParameters, dtype=float)
    index = [0] * nParameters
    def recurse(i):
        if i == nParameters:
            grid[tuple(index)] = viewer.update(parameters)
        else:
            for j, v in enumerate(ranges[i]):
                index[i] = j
                parameters[i] = v
                recurse(i + 1)
    recurse(0)
    return numpy.exp(-grid)

def marginalize(grid, i, j):
    mgrid = grid
    lshape = list(grid.shape)
    for k in range(len(lshape)):
        if k != i and k != j:
            lshape[k] = 1
            mgrid = mgrid.sum(axis=k).reshape(*lshape)
    return mgrid.squeeze()

def plot(grid, viewer, ranges):
    nParameters = len(ranges)
    mid = lambda x: 0.5*(x[:-1] + x[1:])
    for i in range(nParameters):
        for j in range(i + 1):
            pyplot.subplot(nParameters, nParameters, i * nParameters + j + 1)
            mgrid = marginalize(grid, i, j)
            if i == j:
                pyplot.plot(ranges[i], mgrid)
                pyplot.xlim(ranges[i][0], ranges[i][-1])
            else:
                pyplot.contourf(ranges[j], ranges[i], mgrid.transpose())
                pyplot.xlim(ranges[j][0], ranges[j][-1])
                pyplot.ylim(ranges[i][0], ranges[i][-1])
