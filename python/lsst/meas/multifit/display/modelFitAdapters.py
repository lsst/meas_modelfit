#
# LSST Data Management System
# Copyright 2008-2013 LSST Corporation.
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

import numpy

__all__ = ("ModelFitDataAdapter",)

class ModelFitDataAdapter(object):

    def __init__(self, record):
        self.record = record
        self.pdf = record.getPdf()
        self.samples = record.getSamples().copy(deep=True)
        self.values = self.samples["parameters"]
        self.weights = self.samples["weight"]
        self.dimensions = list(record.getInterpreter().getParameterNames())
        self.setRangesFromQuantiles(0.001, 0.999)
        assert self.values.shape[1] == len(self.dimensions)

    def setRangesFromQuantiles(self, lower, upper):
        fractions = numpy.array([lower, upper], dtype=float)
        ranges = self.record.getInterpreter().computeParameterQuantiles(self.record, fractions)
        self.lower = {dim: ranges[i,0] for i, dim in enumerate(self.dimensions)}
        self.upper = {dim: ranges[i,1] for i, dim in enumerate(self.dimensions)}

    def eval1d(self, dim, x):
        i = self.dimensions.index(dim)
        z = numpy.zeros(x.shape, dtype=float)
        if i >= self.pdf.getDimension():
            return None
        projection = self.pdf.project(i)
        projection.evaluate(x.reshape(x.shape + (1,)), z)
        return z

    def eval2d(self, xDim, yDim, x, y):
        i = self.dimensions.index(yDim)
        j = self.dimensions.index(xDim)
        z = numpy.zeros(x.size, dtype=float)
        if i >= self.pdf.getDimension() or j >= self.pdf.getDimension():
            return None
        projection = self.pdf.project(j, i)
        xy = numpy.zeros((x.size, 2), dtype=float)
        xy[:,0] = x.flatten()
        xy[:,1] = y.flatten()
        projection.evaluate(xy, z)
        return z.reshape(x.shape)
