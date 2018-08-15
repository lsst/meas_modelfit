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
from .densityPlot import mergeDefaults
from .. import modelfitLib

__all__ = ("SamplingDataAdapter", "OptimizerTrackLayer", "OptimizerDataAdapter",)


class ModelFitDataAdapter:

    def __init__(self, record):
        self.record = record
        self.pdf = record.getPdf()
        self.dimensions = list(record.getInterpreter().getParameterNames())

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
        xy[:, 0] = x.flatten()
        xy[:, 1] = y.flatten()
        projection.evaluate(xy, z)
        return z.reshape(x.shape)


class SamplingDataAdapter(ModelFitDataAdapter):

    def __init__(self, record):
        ModelFitDataAdapter.__init__(self, record)
        self.samples = record.getSamples().copy(deep=True)
        self.values = self.samples["parameters"]
        self.weights = self.samples["weight"]
        self.setRangesFromQuantiles(0.001, 0.999)
        assert self.values.shape[1] == len(self.dimensions)

    def setRangesFromQuantiles(self, lower, upper):
        fractions = numpy.array([lower, upper], dtype=float)
        ranges = self.record.getInterpreter().computeParameterQuantiles(self.record, fractions)
        self.lower = {dim: ranges[i, 0] for i, dim in enumerate(self.dimensions)}
        self.upper = {dim: ranges[i, 1] for i, dim in enumerate(self.dimensions)}


class OptimizerTrackLayer:

    defaults = dict(
        accepted=dict(
            marker='.', linestyle='-', color='c',
            markevery=(1, 1),  # (start, stride): don't put a marker on the first point
        ),
        rejected=dict(
            marker='.', linestyle='-', color='k', alpha=0.5,
            markevery=3,  # marker at every third point, so we only mark the rejected points
        ),
    )

    def __init__(self, tag, accepted=None, rejected=None):
        self.tag = tag
        self.accepted = mergeDefaults(accepted, self.defaults['accepted'])
        self.rejected = mergeDefaults(rejected, self.defaults['rejected'])

    def plotX(self, axes, data, dim):
        pass

    def plotY(self, axes, data, dim):
        pass

    def plotXY(self, axes, data, xDim, yDim):
        i = data.dimensions.index(yDim)
        j = data.dimensions.index(xDim)
        artists = []
        artists.extend(axes.plot(data.rejected[:, j], data.rejected[:, i], **self.rejected))
        artists.extend(axes.plot(data.accepted[:, j], data.accepted[:, i], **self.accepted))
        return artists


class OptimizerDataAdapter(ModelFitDataAdapter):

    def __init__(self, record):
        ModelFitDataAdapter.__init__(self, record)
        self.samples = record.getSamples().copy(deep=True)
        self.parameters = self.samples["parameters"]
        self.state = self.samples["state"]
        # The first point is neither accepted nor rejected, so we test on rejected and !rejected so
        # as to include the first point with the accepted points
        mask = (self.state & modelfitLib.Optimizer.STATUS_STEP_REJECTED).astype(bool)
        self.accepted = self.parameters[numpy.logical_not(mask)]
        # For each rejected point, we have three path points: the rejected point, the last accepted point,
        # and a NaN to tell matplotlib not to connect to the next one.
        # Note that the defaults for OptimizerTrackLayer use markevery=3 to only put markers on
        # the rejected points
        rejected = []
        current = self.parameters[0]
        nans = numpy.array([numpy.nan] * self.parameters.shape[1], dtype=float)
        for parameters, isRejected in zip(self.parameters, mask):
            if isRejected:
                rejected.extend([parameters, current, nans])
            else:
                current = parameters
        self.rejected = numpy.array(rejected)
        self.lower = {}
        self.upper = {}
        for i, dim in enumerate(self.dimensions):
            projected = self.pdf[0].project(i)
            mu = projected.getMu()
            sigma = projected.getSigma()**0.5
            self.lower[dim] = min(self.accepted[:, i].min(), mu - 3*sigma)
            self.upper[dim] = max(self.accepted[:, i].max(), mu + 3*sigma)
        # Now we setup some special points for a CrossPointsLayer
        self.points = numpy.zeros((2, self.parameters.shape[1]), dtype=float)
        record.getInterpreter().packParameters(
            self.record['initial.nonlinear'], self.record['initial.amplitudes'],
            self.points[0, :]
        )
        record.getInterpreter().packParameters(
            self.record['fit.nonlinear'], self.record['fit.amplitudes'],
            self.points[1, :]
        )
