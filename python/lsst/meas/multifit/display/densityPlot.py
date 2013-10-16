import numpy
import matplotlib

def hide_xticklabels(axes):
    for label in axes.get_xticklabels():
        label.set_visible(False)

def hide_yticklabels(axes):
    for label in axes.get_yticklabels():
        label.set_visible(False)

class TestData(object):

    def __init__(self):
        self.dimensions = ["a", "b", "c"]
        self.ranges = numpy.array([[-9.0, 9.0],
                                   [-6.0, 6.0],
                                   [-3.0, 3.0]], dtype=float)
        self.values = numpy.random.randn(500, 3)
        self.values[:,0] *= 3.0
        self.values[:,1] *= 2.0

    def hist1d(self, dim, limits, bins):
        i = self.dimensions.index(dim)
        return numpy.histogram(self.values[:,i], bins=bins, range=limits, normed=True)

class HistogramLayer(object):

    def __init__(self, bins1d, bins2d, kwds1d=None, kwds2d=None):
        self.bins1d = int(bins1d)
        self.bins2d = int(bins2d)
        self.kwds1d = dict(kwds1d)
        self.kwds2d = dict(kwds2d)

    def hist1d(self, data, dim, limits):
        return data.hist1d(dim, limits, self.bins1d)

    def hist2d(self, data, xDim, yDim, xLimits, yLimits):
        return data.hist2d(xDim, yDim, xLimits, yLimits)

    def plotX(self, axes, data, dim):
        y, xEdge = self.hist1d(data, dim, axes.get_xlim())
        xCenter = 0.5*(xEdge[:-1] + xEdge[1:])
        width = xEdge[1:] - xEdge[:-1]
        axes.bar(xCenter, y, width=width, align='center', **self.kwds1d)

    def plotY(self, axes, data, dim):
        x, yEdge = self.hist1d(data, dim, axes.get_ylim())
        yCenter = 0.5*(yEdge[:-1] + yEdge[1:])
        height = yEdge[1:] - yEdge[:-1]
        axes.barh(yCenter, x, height=height, align='center', **self.kwds1d)

    def plotXY(self, axes, data, xDim, yDim):
        z, xEdge, yEdge = self.hist2d(data, xDim, yDim, axes.get_xlim(), axes.get_ylim())
        axes.imshow(z, aspect='auto', extent=(xEdge[0], xEdge[-1], yEdge[0], yEdge[-1]),
                    origin='lower', **self.kwds2d)

class SurfaceLayer(object):

    def __init__(self, kwds1d, kwds2d):
        self.kwds1d = dict(kwds1d)
        self.kwds2d = dict(kwds2d)

    def eval1d(self, data, dim, x):
        

    def plotX(self, axes, data, dim):
        y, xEdge = self.hist1d(data, dim, axes.get_xlim())
        xCenter = 0.5*(xEdge[:-1] + xEdge[1:])
        width = xEdge[1:] - xEdge[:-1]
        axes.bar(xCenter, y, width=width, align='center', **self.kwds1d)

    def plotY(self, axes, data, dim):
        x, yEdge = self.hist1d(data, dim, axes.get_ylim())
        yCenter = 0.5*(yEdge[:-1] + yEdge[1:])
        height = yEdge[1:] - yEdge[:-1]
        axes.barh(yCenter, x, height=height, align='center', **self.kwds1d)

    def plotXY(self, axes, data, xDim, yDim):
        z, xEdge, yEdge = self.hist2d(data, xDim, yDim, axes.get_xlim(), axes.get_ylim())
        axes.imshow(z, aspect='auto', extent=(xEdge[0], xEdge[-1], yEdge[0], yEdge[-1]),
                    origin='lower', **self.kwds2d)

class DensityPlot(object):

    def __init__(self, figure, data):
        self.figure = figure
        self.dimensions = tuple(data.dimensions)
        self._active = self.dimensions
        self._all_dims = frozenset(self._active)
        if len(self._all_dims) != len(self._active):
            raise ValueError("Dimensions list contains duplicates")
        self.figure.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95, hspace=0.0, wspace=0.0)
        self.ranges = numpy.array(data.ranges)
        self._build_axes()
        self._layers = []

    @property
    def layers(self): return tuple(self._layers)

    def addLayer(self, layer):
        self._layers.append(layer)
        self._replotLayer(layer)

    def _replotLayer(self, layer):
        for i, yDim in self._active:
            layer.plotX(self._axes[None,i], self.data, yDim)
            layer.plotY(self._axes[i,None], self.yDim)
            for j, xDim in self._active:
                if i == j: continue
                layer.plotXY(self._axes[i,j], self.data, xDim, yDim)

    def _replotBox(self, xDim, yDim):
        i = self._active.index(yDim)
        j = self._active.index(xDim)
        for z, layer in enumerate(self._layers):
            layer.plotX(self._axes[None,j], self.data, xDim)
            layer.plotY(self._axes[i,None], self.data, yDim)
            layer.plotXY(self._axes[i,j], self.data, xDim, yDim)

    def _get_active(self): return self._active
    def _set_active(self, active):
        s = set(active)
        if len(s) != len(active):
            raise ValueError("Active set contains duplicates")
        if not self._all_dims.issuperset(s):
            raise ValueError("Invalid values in active set")
        self._active = tuple(active)
    active = property(_get_active, _set_active, doc="sequence of active dimensions to plot")

    def _build_axes(self):
        self._axes = dict()
        n = len(self._active) + 1
        for i in range(n - 1):
            self._axes[None,i] = self.figure.add_subplot(n, n, i+2)
            self._axes[None,i].xaxis.tick_top()
            self._axes[None,i].set_xlim(self.ranges[i,0], self.ranges[i,1])
            hide_yticklabels(self._axes[None,i])
            bbox = self._axes[None,i].get_position()
            bbox.y0 += 0.025
            self._axes[None,i].set_position(bbox)
            self._axes[i,None] = self.figure.add_subplot(n, n, (i+1)*n + 1)
            self._axes[i,None].yaxis.tick_left()
            self._axes[i,None].set_xlim(self.ranges[i,0], self.ranges[i,1])
            hide_xticklabels(self._axes[i,None])
            bbox = self._axes[i,None].get_position()
            bbox.x1 -= 0.025
            self._axes[i,None].set_position(bbox)
        for i in range(n - 1):
            for j in range(i + 1, n - 1):
                self._axes[i,j] = self.figure.add_subplot(n, n, (i+1)*n+j+2,
                                                          sharex=self._axes[None,j],
                                                          sharey=self._axes[i,None])
                self._axes[j,i] = self.figure.add_subplot(n, n, (j+1)*n+i+2,
                                                          sharex=self._axes[None,i],
                                                          sharey=self._axes[j,None])
                hide_yticklabels(self._axes[i,j])
                hide_xticklabels(self._axes[i,j])
                hide_yticklabels(self._axes[j,i])
                hide_xticklabels(self._axes[j,i])

    def draw(self):
        self.figure.canvas.draw()

def replot():
    fig = matplotlib.pyplot.figure()
    return DensityPlot(fig, 3, TestData())
