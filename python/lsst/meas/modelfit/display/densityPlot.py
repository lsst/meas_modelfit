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

"""A set of matplotlib-based classes that displays a grid of 1-d and 2-d slices through an
N-d density.

The main class, DensityPlot, manages the grid of matplotlib.axes.Axes objects, and holds
a sequence of Layer objects that each know how to draw individual 1-d or 2-d plots and a
data object that abstracts away how the N-d density data is actually represented.

For simple cases, users can just create a custom data class with an interface like that of
the ExampleData class provided here, and use the provided HistogramLayer and SurfaceLayer
classes directly.  In more complicated cases, users may want to create their own Layer classes,
which may define their own relationship with the data object.
"""
from builtins import range
from builtins import object

import collections
import numpy
import matplotlib

__all__ = ("HistogramLayer", "SurfaceLayer", "ScatterLayer", "CrossPointsLayer",
           "DensityPlot", "ExampleData", "demo")


def hide_xticklabels(axes):
    for label in axes.get_xticklabels():
        label.set_visible(False)


def hide_yticklabels(axes):
    for label in axes.get_yticklabels():
        label.set_visible(False)


def mergeDefaults(kwds, defaults):
    copy = defaults.copy()
    if kwds is not None:
        copy.update(**kwds)
    return copy


class HistogramLayer(object):
    """A Layer class for DensityPlot for gridded histograms, drawing bar plots in 1-d and
    colormapped large-pixel images in 2-d.

    Relies on two data object attributes:

       values ----- a (M,N) array of data points, where N is the dimension of the dataset and M is the
                    number of data points

       weights ---- (optional) an array of weights with shape (M,); if not present, all weights will
                    be set to unity

    The need for these data object attributes can be removed by subclassing HistogramLayer and overriding
    the hist1d and hist2d methods.
    """

    defaults1d = dict(facecolor='b', alpha=0.5)
    defaults2d = dict(cmap=matplotlib.cm.Blues, vmin=0.0, interpolation='nearest')

    def __init__(self, tag, bins1d=20, bins2d=(20, 20), kwds1d=None, kwds2d=None):
        self.tag = tag
        self.bins1d = bins1d
        self.bins2d = bins2d
        self.kwds1d = mergeDefaults(kwds1d, self.defaults1d)
        self.kwds2d = mergeDefaults(kwds2d, self.defaults2d)

    def hist1d(self, data, dim, limits):
        """Extract points from the data object and compute a 1-d histogram.

        Return value should match that of numpy.histogram: a tuple of (hist, edges),
        where hist is a 1-d array with size=bins1d, and edges is a 1-d array with
        size=self.bins1d+1 giving the upper and lower edges of the bins.
        """
        i = data.dimensions.index(dim)
        if hasattr(data, "weights") and data.weights is not None:
            weights = data.weights
        else:
            weights = None
        return numpy.histogram(data.values[:, i], bins=self.bins1d, weights=weights,
                               range=limits, normed=True)

    def hist2d(self, data, xDim, yDim, xLimits, yLimits):
        """Extract points from the data object and compute a 1-d histogram.

        Return value should match that of numpy.histogram2d: a tuple of (hist, xEdges, yEdges),
        where hist is a 2-d array with shape=bins2d, xEdges is a 1-d array with size=bins2d[0]+1,
        and yEdges is a 1-d array with size=bins2d[1]+1.
        """
        i = data.dimensions.index(yDim)
        j = data.dimensions.index(xDim)
        if hasattr(data, "weights") and data.weights is not None:
            weights = data.weights
        else:
            weights = None
        return numpy.histogram2d(data.values[:, j], data.values[:, i], bins=self.bins2d, weights=weights,
                                 range=(xLimits, yLimits), normed=True)

    def plotX(self, axes, data, dim):
        y, xEdge = self.hist1d(data, dim, axes.get_xlim())
        xCenter = 0.5*(xEdge[:-1] + xEdge[1:])
        width = xEdge[1:] - xEdge[:-1]
        return axes.bar(xCenter, y, width=width, align='center', **self.kwds1d)

    def plotY(self, axes, data, dim):
        x, yEdge = self.hist1d(data, dim, axes.get_ylim())
        yCenter = 0.5*(yEdge[:-1] + yEdge[1:])
        height = yEdge[1:] - yEdge[:-1]
        return axes.barh(yCenter, x, height=height, align='center', **self.kwds1d)

    def plotXY(self, axes, data, xDim, yDim):
        z, xEdge, yEdge = self.hist2d(data, xDim, yDim, axes.get_xlim(), axes.get_ylim())
        return axes.imshow(z.transpose(), aspect='auto', extent=(xEdge[0], xEdge[-1], yEdge[0], yEdge[-1]),
                           origin='lower', **self.kwds2d)


class ScatterLayer(object):
    """A Layer class that plots individual points in 2-d, and does nothing in 1-d.

    Relies on two data object attributes:

       values ----- a (M,N) array of data points, where N is the dimension of the dataset and M is the
                    number of data points

       weights ---- (optional) an array of weights with shape (M,); will be used to set the color of points

    """

    defaults = dict(linewidth=0, alpha=0.2)

    def __init__(self, tag, **kwds):
        self.tag = tag
        self.kwds = mergeDefaults(kwds, self.defaults)

    def plotX(self, axes, data, dim):
        pass

    def plotY(self, axes, data, dim):
        pass

    def plotXY(self, axes, data, xDim, yDim):
        i = data.dimensions.index(yDim)
        j = data.dimensions.index(xDim)
        if hasattr(data, "weights") and data.weights is not None:
            args = data.values[:, j], data.values[:, i], data.weights
        else:
            args = data.values[:, j], data.values[:, i]
        return axes.scatter(*args, **self.kwds)


class SurfaceLayer(object):
    """A Layer class for analytic N-d distributions that can be evaluated in 1-d or 2-d slices.

    The 2-d slices are drawn as contours, and the 1-d slices are drawn as simple curves.

    Relies on eval1d and eval2d methods in the data object; this can be avoided by subclassing
    SurfaceLayer and reimplementing its own eval1d and eval2d methods.
    """

    defaults1d = dict(linewidth=2, color='r')
    defaults2d = dict(linewidths=2, cmap=matplotlib.cm.Reds)

    def __init__(self, tag, steps1d=200, steps2d=200, filled=False, kwds1d=None, kwds2d=None):
        self.tag = tag
        self.steps1d = int(steps1d)
        self.steps2d = int(steps2d)
        self.filled = bool(filled)
        self.kwds1d = mergeDefaults(kwds1d, self.defaults1d)
        self.kwds2d = mergeDefaults(kwds2d, self.defaults2d)

    def eval1d(self, data, dim, x):
        """Return analytic function values for the given values."""
        return data.eval1d(dim, x)

    def eval2d(self, data, xDim, yDim, x, y):
        """Return analytic function values for the given values."""
        return data.eval2d(xDim, yDim, x, y)

    def plotX(self, axes, data, dim):
        xMin, xMax = axes.get_xlim()
        x = numpy.linspace(xMin, xMax, self.steps1d)
        z = self.eval1d(data, dim, x)
        if z is None:
            return
        return axes.plot(x, z, **self.kwds1d)

    def plotY(self, axes, data, dim):
        yMin, yMax = axes.get_ylim()
        y = numpy.linspace(yMin, yMax, self.steps1d)
        z = self.eval1d(data, dim, y)
        if z is None:
            return
        return axes.plot(z, y, **self.kwds1d)

    def plotXY(self, axes, data, xDim, yDim):
        xMin, xMax = axes.get_xlim()
        yMin, yMax = axes.get_ylim()
        xc = numpy.linspace(xMin, xMax, self.steps2d)
        yc = numpy.linspace(yMin, yMax, self.steps2d)
        xg, yg = numpy.meshgrid(xc, yc)
        z = self.eval2d(data, xDim, yDim, xg, yg)
        if z is None:
            return
        if self.filled:
            return axes.contourf(xg, yg, z, 6, **self.kwds2d)
        else:
            return axes.contour(xg, yg, z, 6, **self.kwds2d)


class CrossPointsLayer(object):
    """A layer that marks a few points with axis-length vertical and horizontal lines.

    This relies on a "points" data object attribute.
    """

    defaults = dict(alpha=0.8)

    def __init__(self, tag, colors=("y", "m", "c", "r", "g", "b"), **kwds):
        self.tag = tag
        self.colors = colors
        self.kwds = mergeDefaults(kwds, self.defaults)

    def plotX(self, axes, data, dim):
        i = data.dimensions.index(dim)
        artists = []
        for n, point in enumerate(data.points):
            artists.append(axes.axvline(point[i], color=self.colors[n % len(self.colors)], **self.kwds))
        return artists

    def plotY(self, axes, data, dim):
        i = data.dimensions.index(dim)
        artists = []
        for n, point in enumerate(data.points):
            artists.append(axes.axhline(point[i], color=self.colors[n % len(self.colors)], **self.kwds))
        return artists

    def plotXY(self, axes, data, xDim, yDim):
        i = data.dimensions.index(yDim)
        j = data.dimensions.index(xDim)
        artists = []
        for n, point in enumerate(data.points):
            artists.append(axes.axvline(point[j], color=self.colors[n % len(self.colors)], **self.kwds))
            artists.append(axes.axhline(point[i], color=self.colors[n % len(self.colors)], **self.kwds))
        return artists


class DensityPlot(object):
    """An object that manages a matrix of matplotlib.axes.Axes objects that represent a set of 1-d and 2-d
    slices through an N-d density.
    """

    class LayerDict(collections.MutableMapping):

        def __init__(self, parent):
            self._dict = dict()
            self._parent = parent

        def __delitem__(self, name):
            layer = self._dict.pop(name)
            self._parent._dropLayer(name, layer)

        def __setitem__(self, name, layer):
            self.pop(name, None)
            self._dict[name] = layer
            self._parent._plotLayer(name, layer)

        def __getitem__(self, name):
            return self._dict[name]

        def __iter__(self):
            return iter(self._dict)

        def __len__(self):
            return len(self._dict)

        def __str__(self):
            return str(self._dict)

        def __repr__(self):
            return repr(self._dict)

        def replot(self, name):
            layer = self._dict[name]
            self._parent._dropLayer(name, layer)
            self._parent._plotLayer(name, layer)

    def __init__(self, figure, **kwds):
        self.figure = figure
        self.data = dict(kwds)
        active = []
        self._lower = dict()
        self._upper = dict()
        # We merge the dimension name lists manually rather than using sets to preserve the order.
        # Most of the time we expect all data objects to have the same dimensions anyway.
        for v in self.data.values():
            for dim in v.dimensions:
                if dim not in active:
                    active.append(dim)
                    self._lower[dim] = v.lower[dim]
                    self._upper[dim] = v.upper[dim]
                else:
                    self._lower[dim] = min(v.lower[dim], self._lower[dim])
                    self._upper[dim] = max(v.upper[dim], self._upper[dim])
        self._active = tuple(active)
        self._all_dims = frozenset(self._active)
        self.figure.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95, hspace=0.01, wspace=0.01)
        self._build_axes()
        self.layers = self.LayerDict(self)

    def _dropLayer(self, name, layer):
        def removeArtist(*key):
            try:
                self._objs.pop(key).remove()
            except AttributeError:
                # sometimes the value might be None, which doesn't have a remove
                pass
            except TypeError:
                # probably a matplotlib bug: remove sometimes raises an exception,
                # but it still works
                pass
        for i, yDim in enumerate(self._active):
            removeArtist(None, i, name)
            removeArtist(i, None, name)
            for j, xDim in enumerate(self._active):
                if i == j:
                    continue
                removeArtist(i, j, name)

    def _plotLayer(self, name, layer):
        for i, yDim in enumerate(self._active):
            if yDim not in self.data[layer.tag].dimensions:
                continue
            self._objs[None, i, name] = layer.plotX(self._axes[None, i], self.data[layer.tag], yDim)
            self._objs[i, None, name] = layer.plotY(self._axes[i, None], self.data[layer.tag], yDim)
            for j, xDim in enumerate(self._active):
                if xDim not in self.data[layer.tag].dimensions:
                    continue
                if i == j:
                    continue
                self._objs[i, j, name] = layer.plotXY(self._axes[i, j], self.data[layer.tag], xDim, yDim)
            self._axes[None, i].xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=5, prune='both'))
            self._axes[i, None].yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=5, prune='both'))
            self._axes[None, i].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            self._axes[i, None].yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())

    def _get_active(self): return self._active

    def _set_active(self, active):
        s = set(active)
        if len(s) != len(active):
            raise ValueError("Active set contains duplicates")
        if not self._all_dims.issuperset(s):
            raise ValueError("Invalid values in active set")
        self._active = tuple(active)
        self.replot()
    active = property(_get_active, _set_active, doc="sequence of active dimensions to plot (sequence of str)")

    def replot(self):
        self._lower = {dim: min(self.data[k].lower[dim] for k in self.data) for dim in self._active}
        self._upper = {dim: max(self.data[k].upper[dim] for k in self.data) for dim in self._active}
        self._build_axes()
        for name, layer in self.layers.items():
            self._plotLayer(name, layer)

    def _build_axes(self):
        self.figure.clear()
        self._axes = dict()
        self._objs = dict()
        n = len(self._active)
        iStride = n + 1
        jStride = -1
        iStart = n + 1
        jStart = n
        for i in range(n):
            j = i
            axesX = self._axes[None, j] = self.figure.add_subplot(n+1, n+1, jStart+j*jStride)
            axesX.autoscale(False, axis='x')
            axesX.xaxis.tick_top()
            axesX.set_xlim(self._lower[self._active[j]], self._upper[self._active[j]])
            hide_yticklabels(axesX)
            bbox = axesX.get_position()
            bbox.y1 -= 0.035
            axesX.set_position(bbox)
            axesY = self._axes[i, None] = self.figure.add_subplot(n+1, n+1, iStart + iStart+i*iStride)
            axesY.autoscale(False, axis='y')
            axesY.yaxis.tick_right()
            axesY.set_ylim(self._lower[self._active[i]], self._upper[self._active[i]])
            hide_xticklabels(axesY)
            bbox = axesY.get_position()
            bbox.x1 -= 0.035
            axesY.set_position(bbox)
        for i in range(n):
            for j in range(n):
                axesXY = self._axes[i, j] = self.figure.add_subplot(
                    n+1, n+1, iStart+i*iStride + jStart+j*jStride,
                    sharex=self._axes[None, j],
                    sharey=self._axes[i, None]
                )
                axesXY.autoscale(False)
                if j < n - 1:
                    hide_yticklabels(axesXY)
                if i < n - 1:
                    hide_xticklabels(axesXY)
        for i in range(n):
            j = i
            xbox = self._axes[None, j].get_position()
            ybox = self._axes[i, None].get_position()
            self.figure.text(0.5*(xbox.x0 + xbox.x1), 0.5*(ybox.y0 + ybox.y1), self.active[i],
                             ha='center', va='center', weight='bold')
            self._axes[i, j].get_frame().set_facecolor('none')

    def draw(self):
        self.figure.canvas.draw()


class ExampleData(object):
    """An example data object for DensityPlot, demonstrating the necessarity interface.

    There are two levels of requirements for a data object.  First are the attributes
    required by the DensityPlot object itself; these must be present on every data object:

       dimensions ------ a sequence of strings that provide names for the dimensions

       lower ----------- a dictionary of {dimension-name: lower-bound}

       upper ----------- a dictionary of {dimension-name: upper-bound}

    The second level of requirements are those of the Layer objects provided here.  These
    may be absent if the associated Layer is not used or is subclassed to reimplement the
    Layer method that calls the data object method.  Currently, these include:

       eval1d, eval2d -- methods used by the SurfaceLayer class; see their docs for more info

       values ---------- attribute used by the HistogramLayer and ScatterLayer classes, an array
                         with shape (M,N), where N is the number of dimension and M is the number
                         of data points

       weights --------- optional attribute used by the HistogramLayer and ScatterLayer classes,
                         a 1-d array with size=M that provides weights for each data point
    """

    def __init__(self):
        self.dimensions = ["a", "b", "c"]
        self.mu = numpy.array([-10.0, 0.0, 10.0])
        self.sigma = numpy.array([3.0, 2.0, 1.0])
        self.lower = {dim: -3*self.sigma[i] + self.mu[i] for i, dim in enumerate(self.dimensions)}
        self.upper = {dim: 3*self.sigma[i] + self.mu[i] for i, dim in enumerate(self.dimensions)}
        self.values = numpy.random.randn(2000, 3) * self.sigma[numpy.newaxis, :] + self.mu[numpy.newaxis, :]

    def eval1d(self, dim, x):
        """Evaluate the 1-d analytic function for the given dim at points x (a 1-d numpy array;
        this method must be numpy-vectorized).
        """
        i = self.dimensions.index(dim)
        return numpy.exp(-0.5*((x-self.mu[i])/self.sigma[i])**2) / ((2.0*numpy.pi)**0.5 * self.sigma[i])

    def eval2d(self, xDim, yDim, x, y):
        """Evaluate the 2-d analytic function for the given xDim and yDim at points x,y
        (2-d numpy arrays with the same shape; this method must be numpy-vectorized).
        """
        i = self.dimensions.index(yDim)
        j = self.dimensions.index(xDim)
        return (numpy.exp(-0.5*(((x-self.mu[j])/self.sigma[j])**2 + ((y-self.mu[i])/self.sigma[i])**2))
                / (2.0*numpy.pi * self.sigma[j]*self.sigma[i]))


def demo():
    """Create and return a DensityPlot with example data."""
    fig = matplotlib.pyplot.figure()
    p = DensityPlot(fig, primary=ExampleData())
    p.layers['histogram'] = HistogramLayer('primary')
    p.layers['surface'] = SurfaceLayer('primary')
    p.draw()
    return p
