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
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable

from .densityPlot import mergeDefaults, hide_xticklabels, hide_yticklabels
from .. import multifitLib

__all__ = ("OptimizerDisplay", )

class OptimizerIterationDisplay(object):

    def __init__(self, parent, sample):
        self.parent = parent
        self.sample = sample
        # scale unit grid by trust radius
        self.grid = parent.unitGrid * sample.get(parent.recorder.trust)
        # offset grid to center it on the current parameter point
        self.grid += sample.get(parent.recorder.parameters).reshape((1,)*parent.ndim + (parent.ndim,))
        self._objectiveValues = None
        self._objectiveModel = None
        self.rejected = []

    def __getattr__(self, name):
        # look for keys on the recorder and lookup fields for unknown attributes
        return self.sample.get(getattr(self.parent.recorder, name))

    @property
    def objectiveValues(self):
        if self._objectiveValues is None:
            self._objectiveValues = numpy.zeros(self.grid.shape[:-1], dtype=float)
            self.parent.objective.fillObjectiveValueGrid(self.grid.reshape(-1, self.parent.ndim),
                                                         self._objectiveValues.reshape(-1))
        return self._objectiveValues

    @property
    def objectiveModel(self):
        if self._objectiveModel is None:
            self._objectiveModel = numpy.zeros(self.grid.shape[:-1], dtype=float)
            self.parent.recorder.fillObjectiveModelGrid(self.sample,
                                                        self.grid.reshape(-1, self.parent.ndim),
                                                        self._objectiveModel.reshape(-1))
        return self._objectiveModel

class OptimizerDisplay(object):

    def __init__(self, record, objective, steps=11):
        self.recorder = multifitLib.OptimizerHistoryRecorder(record.getSamples().getSchema())
        # len(dimensions) == N in comments below
        self.dimensions = record.getInterpreter().getParameterNames()
        self.ndim = len(self.dimensions)
        self.track = []
        self.objective = objective
        # This creates a array with shape [steps, ..., steps, N] (a total of N+1 array dimensions):
        # this is an N-dimensional grid, with the last dimension of the grid array giving the coordinates
        # of each grid point.
        # We slice mgrid to generate the basic grid, which is [N, steps, ..., steps]
        mgridArgs = (slice(-1.0, 1.0, steps*1j),) * self.ndim
        # We'll index the result of mgrid with these args to make first dimension last
        transposeArgs = tuple(range(1, self.ndim+1) + [0])
        self.unitGrid = numpy.mgrid[mgridArgs].transpose(transposeArgs).copy()
        current = None
        for sample in record.getSamples():
            if sample.get(self.recorder.state) & multifitLib.Optimizer.STATUS_STEP_REJECTED:
                assert current is not None
                current.rejected.append(sample)
                continue
            current = OptimizerIterationDisplay(self, sample)
            self.track.append(current)

    def plot(self, xDim, yDim, n=0):
        return OptimizerDisplayFigure(self, xDim=xDim, yDim=yDim, n=n)

class OptimizerDisplayFigure(object):

    def __init__(self, parent, xDim, yDim, n=0):
        self.parent = parent
        self.xDim = xDim
        self.yDim = yDim
        self.j = self.parent.dimensions.index(self.xDim)
        self.i = self.parent.dimensions.index(self.yDim)
        self.yKey = self.parent.recorder.parameters[self.i]
        self.xKey = self.parent.recorder.parameters[self.j]
        self.zKey = self.parent.recorder.objective
        # grid slice indices corresponding to the dimensions we're plotting
        self.slice2d = [s//2 for s in self.parent.unitGrid.shape[:-1]]
        self.slice2d[self.i] = slice(None)
        self.slice2d[self.j] = slice(None)
        self.slice2d = tuple(self.slice2d)
        self.sliceX = [s//2 for s in self.parent.unitGrid.shape[:-1]]
        self.sliceX[self.j] = slice(None)
        self.sliceX = tuple(self.sliceX)
        self.sliceY = [s//2 for s in self.parent.unitGrid.shape[:-1]]
        self.sliceY[self.i] = slice(None)
        self.sliceY = tuple(self.sliceY)
        self.track = dict(
            x = numpy.array([iteration.sample.get(self.xKey) for iteration in self.parent.track]),
            y = numpy.array([iteration.sample.get(self.yKey) for iteration in self.parent.track]),
            z = numpy.array([iteration.sample.get(self.zKey) for iteration in self.parent.track]),
            )
        self.n = n
        self.figure = matplotlib.pyplot.figure("%s vs %s" % (xDim, yDim), figsize=(16, 8))
        self.figure.subplots_adjust(left=0.025, right=0.975, bottom=0.08, top=0.95, wspace=0.12)
        self.axes3d = self.figure.add_subplot(1, 2, 1, projection='3d')
        self.axes3d.autoscale(False)
        self.axes3d.set_xlabel(self.xDim)
        self.axes3d.set_ylabel(self.yDim)
        self.axes2d = self.figure.add_subplot(1, 2, 2)
        self.axes2d.set_xlabel(self.xDim)
        self.axes2d.set_ylabel(self.yDim)
        self.axes2d.autoscale(False)
        divider = make_axes_locatable(self.axes2d)
        self.axesX = divider.append_axes("top", 1.5, pad=0.1, sharex=self.axes2d)
        self.axesX.autoscale(False, axis='x')
        hide_xticklabels(self.axesX)
        self.axesY = divider.append_axes("right", 1.5, pad=0.1, sharey=self.axes2d)
        self.axesY.autoscale(False, axis='y')
        hide_yticklabels(self.axesY)
        self.artists = []
        self.guessExtent()
        self.plotTrack()
        self.plotRejected()
        self.plotSurfaces()

    @property
    def xlim(self): return self._extent[:2]

    @property
    def ylim(self): return self._extent[2:4]

    @property
    def zlim(self): return self._extent[4:]

    def guessExtent(self):
        current = self.parent.track[self.n]
        x = current.sample.get(self.xKey)
        y = current.sample.get(self.yKey)
        zMin1 = current.objectiveValues[self.slice2d].min()
        zMax1 = current.objectiveValues[self.slice2d].max()
        zMin2 = current.objectiveModel[self.slice2d].min()
        zMax2 = current.objectiveModel[self.slice2d].max()
        self.setExtent(x0=x - current.trust, x1=x + current.trust,
                       y0=y - current.trust, y1=y + current.trust,
                       z0=min(zMin1, zMin2), z1=max(zMax1, zMax2), lock=False)

    def setExtent(self, x0=None, x1=None, y0=None, y1=None, z0=None, z1=None, lock=True):
        if x0 is None: x0 = self._extent[0]
        if x1 is None: x1 = self._extent[1]
        if y0 is None: y0 = self._extent[2]
        if y1 is None: y1 = self._extent[3]
        if z0 is None: z0 = self._extent[4]
        if z1 is None: z1 = self._extent[5]
        self._extent = (x0, x1, y0, y1, z0, z1)
        self._lock = lock
        self.axes3d.set_xlim(*self.xlim)
        self.axes3d.set_ylim(*self.ylim)
        self.axes3d.set_zlim(*self.zlim)
        self.axes2d.set_xlim(*self.xlim)
        self.axes2d.set_ylim(*self.ylim)
        self.axesX.set_ylim(*self.zlim)
        self.axesY.set_xlim(*self.zlim)

    def _clipZ(self, x, y, z):
        # clipping is currently disabled; more trouble than it's worth
        if False:
            mask = numpy.logical_or.reduce([x < self.xlim[0], x > self.xlim[1],
                                            y < self.ylim[0], y > self.ylim[1],
                                            z < self.zlim[0], z > self.zlim[1]],
                                           axis=0)

            z[mask] = numpy.nan
            return numpy.logical_not(mask).astype(int).sum() > 4
        return True

    def _contour(self, axes, *args, **kwds):
        self.artists.extend(axes.contour(*args, **kwds).collections)

    def plotTrack(self):
        kwds = dict(markeredgewidth=0, markerfacecolor='g', color='g', marker='o')
        self.axes3d.plot(self.track['x'], self.track['y'], self.track['z'], **kwds)
        self.axes2d.plot(self.track['x'], self.track['y'], **kwds)
        self.axesX.plot(self.track['x'], self.track['z'], **kwds)
        self.axesY.plot(self.track['z'], self.track['y'], **kwds)

    def plotRejected(self):
        kwds = dict(markeredgewidth=0, markerfacecolor='r', color='r', marker='v')
        current = self.parent.track[self.n]
        cx = current.sample.get(self.xKey)
        cy = current.sample.get(self.yKey)
        cz = current.sample.get(self.zKey)
        for r in current.rejected:
            x = [cx, r.get(self.xKey)]
            y = [cy, r.get(self.yKey)]
            z = [cz, r.get(self.zKey)]
            self.artists.extend(self.axes3d.plot(x, y, z, **kwds))
            self.artists.extend(self.axes2d.plot(x, y, **kwds))
            self.artists.extend(self.axesX.plot(x, z, **kwds))
            self.artists.extend(self.axesY.plot(z, y, **kwds))

    def plotSurfaces(self):
        current = self.parent.track[self.n]

        # Start with 2-d and 3-d surfaces
        x = current.grid[self.slice2d + (self.j,)]
        y = current.grid[self.slice2d + (self.i,)]
        z1 = current.objectiveValues[self.slice2d].copy()
        z2 = current.objectiveModel[self.slice2d].copy()
        norm = matplotlib.colors.Normalize(vmin=self.zlim[0], vmax=self.zlim[1])

        self._contour(self.axes2d, x, y, z1, cmap=matplotlib.cm.spring, norm=norm)
        self._contour(self.axes2d, x, y, z2, cmap=matplotlib.cm.winter, norm=norm)

        # matplotlib doesn't do clipping in 3d, so we'll do that manually
        if self._clipZ(x, y, z1):
            self._contour(self.axes3d, x, y, z1, cmap=matplotlib.cm.spring, norm=norm)
            self.artists.append(self.axes3d.plot_surface(x, y, z1, rstride=1, cstride=1,
                                                         cmap=matplotlib.cm.spring, norm=norm,
                                                         linewidth=0, antialiased=1, alpha=0.5))
        if self._clipZ(x, y, z2):
            self._contour(self.axes3d, x, y, z2, cmap=matplotlib.cm.winter, norm=norm)
            self.artists.append(self.axes3d.plot_surface(x, y, z2, rstride=1, cstride=1,
                                                         cmap=matplotlib.cm.winter, norm=norm,
                                                         linewidth=0, antialiased=1, alpha=0.5))

        # Now the 1-d surfaces
        self.artists.extend(self.axesX.plot(current.grid[self.sliceX + (self.j,)],
                                            current.objectiveValues[self.sliceX], 'm-'))
        self.artists.extend(self.axesX.plot(current.grid[self.sliceX + (self.j,)],
                                            current.objectiveModel[self.sliceX], 'c-'))
        self.artists.extend(self.axesY.plot(current.objectiveValues[self.sliceY],
                                            current.grid[self.sliceY + (self.i,)], 'm-'))
        self.artists.extend(self.axesY.plot(current.objectiveModel[self.sliceY],
                                            current.grid[self.sliceY + (self.i,)], 'c-'))

    def move(self, n):
        self.n = n
        if not self._lock:
            self.guessExtent()
        for artist in self.artists:
            try:
                artist.remove()
            except TypeError:
                # sometimes matplotlib throws an exception even though everything worked fine
                pass
        self.artists = []
        self.plotSurfaces()
        self.plotRejected()
        self.figure.canvas.draw()
