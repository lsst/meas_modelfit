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
from .. import multifitLib

__all__ = ("OptimizerDisplay",)

class OptimizerIterationDisplay(object):

    def __init__(self, parent, sample, ):
        self.parent = parent
        self.sample = sample
        # idx: tuple index equivalent to [:, numpy.newaxis, numpy.newaxis, ...] for all dimensions
        idx = (slice(None),) + (numpy.newaxis,)*len(self.parent.dimensions)
        self.grid = parent.grid + self.sample.get(parent.recorder.parameters)[idx]
        self._objectiveValues = None
        self._objectiveModel = None
        self.rejected = []

    @property
    def objectiveValues(self):
        if self._objectiveValues is None:
            parameterDim = len(self.parent.dimensions)
            self._objectiveValues = numpy.zeros(self.grid.shape[:-1], dtype=float)
            self.parent.objective.fillObjectiveValueGrid(self.grid.reshape(-1, parameterDim),
                                                         self._objectiveValues.reshape(-1))
        return self._objectiveValues

    @property
    def objectiveModel(self):
        if self._objectiveModel is None:
            parameterDim = len(self.parent.dimensions)
            self._objectiveModel = numpy.zeros(self.grid.shape[:-1], dtype=float)
            self.parent.objective.fillObjectiveValueGrid(self.grid.reshape(-1, parameterDim),
                                                         self._objectiveModel.reshape(-1))
        return self._objectiveModel

class OptimizerDisplay(object):

    def __init__(self, record, objective, bound=1.0, steps=21):
        self.recorder = multifitLib.OptimizerHistoryRecorder(record.getSchema())
        # len(dimensions) == N in comments below
        self.dimensions = record.getInterpreter().getParameterNames()
        self.track = []
        self.objective = objective
        # This creates a array with shape [steps, steps, steps, ..., N] (a total of N+1 array dimensions):
        # this is an N-dimensional grid, with the last dimension of the grid array giving the coordinates
        # of each grid point.
        self.grid = numpy.mgrid[(slice(-bound, bound, steps*1j),) * len(self.dimensions)].transpose().copy()
        current = None
        for sample in record.getSamples():
            if sample.get(recorder.state) & multifitLib.Optimizer.STATE_STEP_REJECTED:
                assert current is not None
                current.rejected.append(sample)
                continue
            current = OptimizerIterationDisplay(self, sample)
            self.track.append(current)
            

        # The first point is neither accepted nor rejected, so we test on rejected and !rejected so
        # as to include the first point with main track.
        mask = (self.state & multifitLib.Optimizer.STATUS_STEP_REJECTED).astype(bool)

        self.track = self.parameters[numpy.logical_not(mask)]
