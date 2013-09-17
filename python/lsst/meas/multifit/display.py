#!/usr/bin/env python
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

import os
import collections
import numpy
import copy
from matplotlib import pyplot
import matplotlib.colors
import mpl_toolkits.mplot3d

import lsst.pipe.base
import lsst.pex.config
import lsst.daf.persistence
import lsst.meas.extensions.multiShapelet
import lsst.afw.image
import lsst.afw.geom.ellipses
import lsst.afw.display.ds9

from .measureCcd import *
from .measureCoadd import *
from .measureMulti import *
from .multifitLib import *

__all__ = ("makeHistogramGrid", "InteractiveFitter")

class DensityPlotSet(object):

    plabels = dict(
        Quadrupole=("Ixx", "Iyy", "Ixy"),
        Axes=("a", "b", "theta"),
        )
    elabels = dict(
        ConformalShear="$\\eta",
        ReducedShear="$g",
        Distortion="$e",
        )
    rlabels = dict(
        DeterminantRadius="$r_d$",
        TraceRadius="$r_t$",
        LogDeterminantRadius="$\\ln{r_d}$",
        LogTraceRadius="$\\ln{r_t}$",
        )
    for ek, ev in elabels.items():
        for rk, rv in rlabels.items():
            plabels["Separable%s%s" % (ek, rk)] = (ev + "_1$", ev + "_2$", rv)

    defaults = dict(
        hist1d=dict(facecolor='b'),
        hist2d=dict(cmap=pyplot.cm.Blues, vmin=0.0, interpolation='nearest'),
        scatter3d=dict(linewidth=0, alpha=0.5, cmap=pyplot.cm.Greys, s=4.0,
                       norm=matplotlib.colors.LogNorm(clip=True)),
        scatter2d=dict(linewidth=0, alpha=0.5, cmap=pyplot.cm.Greys, s=4.0,
                       norm=matplotlib.colors.LogNorm(clip=True)),
        mixture1d=dict(linewidth=2, color='r'),
        mixture2d=dict(linewidths=2, cmap=pyplot.cm.Reds, vmin=0.0, levels=6),
        mixture3d=dict(linewidths=2, cmap=pyplot.cm.Reds, vmin=0.0, levels=6),
        )

    def __init__(self, record, samples=None, mixture=None, bins1d=32, bins2d=32, fraction=0.9999, label=None,
                 **kwds):
        """Plot the samples from a ModelFitRecord in histogram and scatter form, possibly with the
        proposal distribution for comparison.  The samples must have three nonlinear dimensions.

        @param[in]  record       ModelFitRecord from which to obtain data (samples, truth values)
        @param[in]  samples      If not None, a SampleSet to use instead of record.getSamples()
        @param[in]  mixture      If not None, a Mixture3 to use instead of record.getProposal()
        @param[in]  bins1d       Number of bins for 1-d histograms.  May be a three-element sequence
                                 to specify different numbers of bins for different dimensions.
        @param[in]  bins2d       Number of bins for 2-d histograms.  May be a three-element sequence
                                 to specify different numbers of bins for different dimensions.
        @param[in]  fraction     Set the bounds such that this (weighted) fraction of the SampleSet
                                 is kept in each dimension, clipping upper and lower evenly.  If bounds
                                 is 0.98, for instance, we'll compute the 1% and 99% percentiles in
                                 the marginal distribution of each dimension separately.
        @param[in]  **kwds       Additional keyword arguments to pass to matplotlib plotting routines,
                                 in the form of nested dictionaries.  For example, passing
                                 "scatter3d=dict(cmap=pyplot.cm.jet)" will set additional keywords
                                 to be passed when making the 3-d scatter plot.  Valid top-level keywords
                                 are scatter2d, scatter3d, hist1d, hist2d, mixture1d, mixture2d, and
                                 mixture3d.  Setting any of these top-level keywords to None will disable
                                 that plot entirely.
        """
        # setup arguments, set defaults
        if samples is None:
            samples = record.getSamples()
        self.samples = samples
        if mixture is None:
            mixture = samples.getProposal()
        if not isinstance(bins2d, collections.Iterable):
            bins2d = [bins2d] * 3
        if not isinstance(bins1d, collections.Iterable):
            bins1d = [bins1d] * 3
        self.kwds = lsst.pipe.base.Struct(**copy.deepcopy(self.defaults))
        for k, v in kwds.iteritems():
            if not hasattr(self.kwds, k):
                raise TypeError("Unrecognized keyword argument: %s" % k)
            if v is None:
                setattr(self.kwds, k, None)
            else:
                d = getattr(self.kwds, k)
                d.update(**v)
        # extract data
        paramDef = samples.getParameterDefinition()
        self.labels = self.plabels[paramDef.name]
        catalog = samples.getCatalog().copy(deep=True)
        self.parameters = catalog["parameters"]
        self.weights = catalog["weight"]
        # extract special values
        self.special = {}
        self.special["mean"] = (samples.computeMean(), "g<")
        refEllipse = paramDef.makeEllipse(self.special["mean"][0]).getCore()
        refEllipse.assign(record.get("ref.ellipse"))
        self.special["truth"] = (refEllipse.getParameterVector(), "yo")
        # setup parameter bounds
        self.ranges = samples.computeQuantiles(numpy.array([0.5-0.5*fraction, 0.5+0.5*fraction])).transpose()
        self.limits = self.ranges.copy()
        self.masks = numpy.ones((3, self.weights.size), dtype=bool)
        for k in range(3):
            self.masks[k] = numpy.logical_and(self.parameters[:,k] >= self.ranges[k,0],
                                              self.parameters[:,k] <= self.ranges[k,1])
        # Compute histograms
        self.histograms = {}
        self.k = [(2,1), (0,2), (0,1)]
        for i, j in self.k:
            h, xe, ye = numpy.histogram2d(self.parameters[:,i], self.parameters[:,j], weights=self.weights,
                                          normed=True, bins=(bins2d[i], bins2d[j]),
                                          range=(self.ranges[i], self.ranges[j]))
            xg, yg = numpy.meshgrid(0.5*(xe[1:]+xe[:-1]), 0.5*(ye[1:]+ye[:-1]))
            self.histograms[i,j] = (h.transpose(), xg, yg)
        for i in range(3):
            h, xe = numpy.histogram(self.parameters[:,i], weights=self.weights,
                                    normed=True, bins=bins1d[i], range=self.ranges[i])
            xg = 0.5*(xe[1:]+xe[:-1])
            w = xe[1:] - xe[:-1]
            self.histograms[i] = (h, xg, w)
        # Evaluate mixture distribution
        self.mixtures = None
        if mixture is not None:
            mixture = Mixture3.cast(mixture)
            self.mixtures = {}
            for i in range(3):
                mix1d = mixture.project(i)
                x = numpy.linspace(self.ranges[i,0], self.ranges[i,1], 100)
                p = numpy.zeros(x.shape, dtype=float)
                c = numpy.zeros((x.shape[0], len(mix1d)), dtype=float)
                mix1d.evaluate(x.reshape(-1,1), p)
                mix1d.evaluateComponents(x.reshape(-1,1), c)
                self.mixtures[i] = (p, c, x)
            for i, j in self.k:
                mix2d = mixture.project(i, j)
                x, y = numpy.meshgrid(self.mixtures[i][2], self.mixtures[j][2])
                xy = numpy.c_[x.flatten(), y.flatten()].copy()
                p = numpy.zeros(xy.shape[0], dtype=float)
                mix2d.evaluate(xy, p)
                self.mixtures[i,j] = (p.reshape(x.shape), x, y)
        # prepare plot axes
        if label is None:
            label = "%s" % record.getId()
        self.figure = pyplot.figure(label, figsize=(9.33, 5.67))
        self.figure.clear()
        self.axes = {}
        self.axes[0,1,2] = self.figure.add_subplot(221, projection='3d')
        self.axes[0,1,2].set_xlabel("\n"+self.labels[0])
        self.axes[0,1,2].set_ylabel("\n"+self.labels[1])
        self.axes[0,1,2].set_zlabel("\n"+self.labels[2])
        self.axes[0,1,2].patch.set_visible(False)
        self.axes[0,1,2].xaxis.set_pane_color((1,1,1))
        self.axes[0,1,2].yaxis.set_pane_color((1,1,1))
        self.axes[0,1,2].zaxis.set_pane_color((1,1,1))
        self.axes[0,1] = self.figure.add_subplot(224)
        self.axes[0,1].yaxis.set_label_position("right")
        self.axes[0,1].yaxis.tick_right()
        self.axes[0,1].xaxis.set_label_position("bottom")
        self.axes[0,1].xaxis.tick_bottom()
        self.axes[0,1].callbacks.connect("xlim_changed", lambda ax: self.setRange(0, *ax.get_xlim()))
        self.axes[0,1].callbacks.connect("ylim_changed", lambda ax: self.setRange(1, *ax.get_ylim()))
        self.axes[0,1].set_xlabel(self.labels[0])
        self.axes[0,1].set_ylabel(self.labels[1])
        self.axes[0,2] = self.figure.add_subplot(222)
        self.axes[0,2].xaxis.set_label_position("top")
        self.axes[0,2].xaxis.tick_top()
        self.axes[0,2].set_xlabel(self.labels[0])
        self.axes[0,2].yaxis.set_label_position("right")
        self.axes[0,2].yaxis.tick_right()
        self.axes[0,2].set_ylabel(self.labels[2])
        self.axes[0,2].callbacks.connect("xlim_changed", lambda ax: self.setRange(0, *ax.get_xlim()))
        self.axes[0,2].callbacks.connect("ylim_changed", lambda ax: self.setRange(2, *ax.get_ylim()))
        self.axes[2,1] = self.figure.add_subplot(223)
        self.axes[2,1].yaxis.set_label_position("left")
        self.axes[2,1].yaxis.tick_left()
        self.axes[2,1].set_ylabel(self.labels[1])
        self.axes[2,1].xaxis.set_label_position("bottom")
        self.axes[2,1].xaxis.tick_bottom()
        self.axes[2,1].set_xlabel(self.labels[2])
        self.axes[2,1].callbacks.connect("xlim_changed", lambda ax: self.setRange(2, *ax.get_xlim()))
        self.axes[2,1].callbacks.connect("ylim_changed", lambda ax: self.setRange(1, *ax.get_ylim()))
        bbox1 = (0.08, 0.08, 0.6, 0.92)
        space1 = 0.05
        self.figure.subplots_adjust(left=bbox1[0], bottom=bbox1[1], right=bbox1[2], top=bbox1[3],
                                    wspace=space1, hspace=space1)
        space2 = space1*1.5
        height2 = ((bbox1[3] - bbox1[1]) - space2*2) / 3
        step2 = height2 + space2
        onUpdate = [lambda ax: self.setRange(0, *ax.get_xlim()),
                    lambda ax: self.setRange(1, *ax.get_xlim()),
                    lambda ax: self.setRange(2, *ax.get_xlim())]
        def replaceTicks(getter, setter):
            ticks = getter()
            labels = [""] * len(ticks)
            labels[0] = ticks[0]
            labels[-1] = ticks[-1]
            setter(labels)
        for i in range(3):
            rect = bbox1[2] + 0.08, bbox1[1]+(2-i)*step2, 0.25, height2
            self.axes[i] = self.figure.add_axes(rect)
            self.axes[i].callbacks.connect("xlim_changed", onUpdate[i])
            self.axes[i].yaxis.tick_right()
            self.axes[i].set_xlabel(self.labels[i])
            replaceTicks(self.axes[0,1,2].get_xticks, self.axes[0,1,2].set_xticklabels)
            replaceTicks(self.axes[0,1,2].get_yticks, self.axes[0,1,2].set_yticklabels)
            replaceTicks(self.axes[0,1,2].get_zticks, self.axes[0,1,2].set_zticklabels)
        for axes in self.axes.itervalues():
            axes.tick_params(labelsize=9)
        # Make plots
        self.update()

    def save(self, filename, **kwds):
        kwds = kwds.copy()
        kwds.setdefault("facecolor", self.figure.get_facecolor())
        self.figure.patch.set_alpha(0.0)
        self.figure.patch.set_facecolor('k') # weirdly needed to save as transparent
        self.figure.savefig(filename, **kwds)

    def update(self):
        for i, j in self.k:
            self.plot2d(i, j)
        for i in range(3):
            self.plot1d(i)
        self.plot3d()
        self.figure.canvas.draw()

    def setRange(self, k, vmin=None, vmax=None, fraction=None, update=True):
        if fraction is not None:
            qr = numpy.array([[0.5-0.5*fraction, 0.5+0.5*fraction]])
            vmin, vmax = self.samples.computeQuantiles(qr, k)
        else:
            if vmin is None: vmin = self.limits[k][0]
            if vmax is None: vmax = self.limits[k][1]
        self.masks[k] = numpy.logical_and(self.parameters[:,k] >= vmin, self.parameters[:,k] <= vmax)
        if k == 0:
            self.axes[0,1].set_xlim(vmin, vmax, emit=False)
            self.axes[0,2].set_xlim(vmin, vmax, emit=False)
        elif k == 1:
            self.axes[0,1].set_ylim(vmin, vmax, emit=False)
            self.axes[2,1].set_ylim(vmin, vmax, emit=False)
        elif k == 2:
            self.axes[0,2].set_ylim(vmin, vmax, emit=False)
            self.axes[2,1].set_xlim(vmin, vmax, emit=False)
        self.axes[k].set_xlim(vmin, vmax, emit=False)
        self.limits[k] = (vmin, vmax)
        if update:
            self.plot3d()
            self.figure.canvas.draw()

    def setContourLevels(self, i, j, kwds):
        if not isinstance(kwds['levels'], collections.Iterable):
            kwds = kwds.copy()
            kwds['levels'] = numpy.linspace(0.0, self.mixtures[i,j][0].max(), kwds['levels']+2)[1:-1]
        return kwds

    def plot1d(self, i):
        del self.axes[i].lines[:]
        del self.axes[i].patches[:]
        if self.kwds.hist1d is not None:
            h, xg, w = self.histograms[i]
            self.axes[i].bar(xg, h, width=w, align='center', **self.kwds.hist1d)
        if self.kwds.mixture1d is not None and self.mixtures is not None:
            h, c, xg = self.mixtures[i]
            self.axes[i].plot(xg, h, **self.kwds.mixture1d)
            for j in range(c.shape[1]):
                self.axes[i].plot(xg, c[:,j], alpha=0.5, **self.kwds.mixture1d)
        for name, (data, symbol) in self.special.items():
            self.axes[i].axvline(data[i], color=symbol[0], linewidth=2)

    def plot2d(self, i, j):
        del self.axes[i,j].lines[:]
        del self.axes[i,j].collections[:]
        del self.axes[i,j].images[:]
        if self.kwds.hist2d is not None:
            h, xg, yg = self.histograms[i,j]
            self.axes[i,j].imshow(h, aspect='auto', extent=(tuple(self.ranges[i])+tuple(self.ranges[j])),
                                  origin='lower', **self.kwds.hist2d)
        if self.kwds.mixture2d is not None and self.mixtures is not None:
            h, xg, yg = self.mixtures[i,j]
            kwds = self.setContourLevels(i, j, self.kwds.mixture2d)
            self.axes[i,j].contour(xg, yg, h, **kwds)
        if self.kwds.scatter2d is not None:
            self.axes[i,j].scatter(self.parameters[:,i], self.parameters[:,j],
                                   c=self.weights, **self.kwds.scatter2d)
        for name, (data, symbol) in self.special.items():
            self.axes[i,j].plot(data[i:i+1], data[j:j+1], symbol, label=name,
                                markersize=9, alpha=0.8)
            self.axes[i,j].axvline(data[i], color=symbol[0])
            self.axes[i,j].axhline(data[j], color=symbol[0])

    def plot3d(self):
        del self.axes[0,1,2].lines[:]
        del self.axes[0,1,2].collections[:]
        if self.kwds.scatter3d is not None:
            mask = numpy.logical_and.reduce(self.masks, axis=0)
            self.axes[0,1,2].scatter(self.parameters[mask,0], self.parameters[mask,1],
                                     self.parameters[mask,2], c=self.weights[mask],
                                     **self.kwds.scatter3d)
        if self.kwds.mixture3d is not None and self.mixtures is not None:
            slices = []
            for k in range(3):
                x = self.mixtures[k][2]
                imin = numpy.searchsorted(x, self.limits[k][0], side='left')
                imax = numpy.searchsorted(x, self.limits[k][1], side='right')
                slices.append(slice(imin,imax+1))
            for k, (i,j) in enumerate(self.k):
                h, xg, yg = self.mixtures[i,j]
                zdir = ['x', 'y', 'z']
                args = [None, None, None]
                args[i] = xg[slices[j], slices[i]]
                args[j] = yg[slices[j], slices[i]]
                args[k] = h[slices[j], slices[i]]
                kwds = self.setContourLevels(i, j, self.kwds.mixture3d)
                offset = self.limits[k][1] if k == 1 else self.limits[k][0]
                self.axes[0,1,2].contour(args[0], args[1], args[2], zdir=zdir[k],
                                         offset=offset, **kwds)
        for name, (data, symbol) in self.special.items():
            self.axes[0,1,2].plot(data[0:1], data[1:2], data[2:3], symbol, label=name,
                                  markersize=9, alpha=0.8)
            self.axes[0,1,2].plot([data[0]]*2, [data[1]]*2, self.limits[2], symbol[0])
            self.axes[0,1,2].plot([data[0]]*2, self.limits[1], [data[2]]*2, symbol[0])
            self.axes[0,1,2].plot(self.limits[0], [data[1]]*2, [data[2]]*2, symbol[0])
        self.axes[0,1,2].set_xlim(*self.limits[0])
        self.axes[0,1,2].set_ylim(*self.limits[1])
        self.axes[0,1,2].set_zlim(*self.limits[2])

class Interactive(object):
    """Interactive analysis helper class

    This class manages a butler, calexp, modelfits catalog, and an instance
    of a Measure*Task, allowing individual objects to be re-fit and plotted.
    """

    def __init__(self, rerun, config=None, dataId=None, mode="ccd"):
        """Construct an interactive analysis object.

        @param[in]  rerun    Output directory, relative to $S13_DATA_DIR/output.
                             measureCcd.py (or measureCoadd.py if coadd=True) must
                             have been run (possibly with prepOnly=True) previously
                             with this output directory.
        @param[in]  config   MeasureCcdTask.ConfigClass instance; if None, it
                             will be loaded from disk.
        @param[in]  dataId   Butler data ID of the image to analyze.
        @param[in]  mode     One of "ccd", "coadd", or "multi", indicating whether
                             to use MeasureCcdTask, MeasureCoaddTask, or MeasureMultiTask.
        """
        root = os.environ["S13_DATA_DIR"]
        self.butler = lsst.daf.persistence.Butler(os.path.join(root, "output", rerun))
        self.mode = mode.lower()
        if self.mode == "ccd":
            if dataId is None:
                dataId = dict(visit=100, raft="2,2", sensor="1,1")
            self.dataRef = self.butler.dataRef("calexp", dataId=dataId)
            if config is None:
                config = self.butler.get("measureCcd_config", immediate=True)
            self.exposure = self.dataRef.get("calexp", immediate=True)
            self.modelfits = self.dataRef.get("modelfits", immediate=True)
            self.task = MeasureCcdTask(config=config)
        elif self.mode == "coadd":
            if dataId is None:
                dataId = dict(tract=0, patch="2,2", filter="i")
            self.dataRef = self.butler.dataRef("deepCoadd_calexp", dataId=dataId)
            if config is None:
                config = self.butler.get("deep_measureCoadd_config", immediate=True)
            self.exposure = self.dataRef.get("deepCoadd_calexp", immediate=True)
            self.modelfits = self.dataRef.get("deepCoadd_modelfits", immediate=True)
            self.task = MeasureCoaddTask(config=config)
        elif self.mode.startswith("multi"):
            if dataId is None:
                dataId = dict(tract=0, patch="2,2", filter="i")
            self.dataRef = self.butler.dataRef("deepCoadd_calexp", dataId=dataId)
            if config is None:
                config = self.butler.get("deep_measureMulti_config", immediate=True)
            self.task = MeasureMultiTask(config=config)
            inputs = self.task.readInputs(self.dataRef)
            self.exposure = inputs.coadd
            self.coaddInputCat = inputs.coaddInputCat
            self.modelfits = self.dataRef.get("deepCoadd_multiModelfits", immediate=True)

    def fit(self, index=0, id=None, doWarmStart=None):
        """Re-fit the object indicated by the given record sequential index
        or source ID, returning the Struct object returned by
        MeasureImageTask.processObject.
        """
        if doWarmStart is None:
            doWarmStart = True if self.mode.startswith('multi') else False
        if id is not None:
            record = self.modelfits.find(id)
        else:
            record = self.modelfits[index]
        if self.mode.startswith("multi"):
            return self.task.processObject(record=record, coadd=self.exposure,
                                           coaddInputCat=self.coaddInputCat)
        else:
            return self.task.processObject(self.exposure, record, doWarmStart=doWarmStart)

    def plotSamples(self, r, iteration=None, suffix="", **kwds):
        """Plot the sammples and proposal distribution from an interactive fit.
        """
        if isinstance(r, ModelFitRecord):
            record = r
            label = "%s%s" % (record.getId(), suffix)
            samples = record.getSamples()
            iterations = []
        else:
            record = r.record
            if iteration is not None:
                samples = r.sampler.iterations[iteration]
                label = "%s+%s%s" % (record.getId(), iteration, suffix)
            else:
                label = "%s%s" % (record.getId(), suffix)
                samples = record.getSamples()
            iterations = r.sampler.iterations
        print "iterations=%s, perplexity=%f, essf=%f" % (len(iterations),
                                                         samples.computeNormalizedPerplexity(),
                                                         samples.computeEffectiveSampleSizeFraction())
        return DensityPlotSet(record, samples=samples, label=label, **kwds)

    def displayResiduals(self, r, parameters=None):
        """Display the data postage stamp along with the model image and residuals in ds9.

        @param[in] r            Result Struct returned by fit()
        @param[in] parameters   Parameter vector for the model; defaults to the posterior mean.
        """
        center = r.record.getPointD(self.task.keys["source.center"])
        if parameters == "truth":
            key = self.task.keys["ref.ellipse"]
            ellipse = lsst.afw.geom.ellipses.Ellipse(r.record.get(key), center)
        else:
            samples = r.record.getSamples()
            if parameters == "mean" or parameters is None:
                parameters = samples.computeMean()
            elif parameters == "median":
                parameters = samples.computeQuantiles(numpy.array([0.5]))
            else:
                parameters = numpy.array(parameters)
            paramDef = samples.getParameterDefinition()
            ellipse = paramDef.makeEllipse(parameters, center)
        joint = r.objective.evaluate(ellipse)
        bbox = r.record.getFootprint().getBBox()
        bbox.grow(2)
        xe = numpy.arange(bbox.getBeginX(), bbox.getEndX(), dtype=numpy.float32)
        ye = numpy.arange(bbox.getBeginY(), bbox.getEndY(), dtype=numpy.float32)
        xg, yg = numpy.meshgrid(xe, ye)
        matrixBuilder = lsst.shapelet.MultiShapeletMatrixBuilderF(
            self.task.basis, r.psf, xg.ravel(), yg.ravel(), self.task.config.objective.useApproximateExp
            )
        matrix = numpy.zeros((joint.grad.size, xg.size), dtype=numpy.float32).transpose()
        matrixBuilder.build(matrix, ellipse)
        mu,_,_,_ = numpy.linalg.lstsq(joint.fisher, -joint.grad)
        array = numpy.dot(matrix, mu).reshape(*xg.shape)

        data = lsst.afw.image.MaskedImageF(self.exposure.getMaskedImage(), bbox,
                                           lsst.afw.image.PARENT, True)
        bitmask = data.getMask().addMaskPlane("FIT_REGION")
        regionMask = lsst.afw.image.MaskU(bbox)
        lsst.afw.detection.setMaskFromFootprint(regionMask, r.record.getFootprint(), bitmask)
        dataMask = data.getMask()
        dataMask |= regionMask
        model = lsst.afw.image.MaskedImageF(lsst.afw.image.ImageF(bbox), regionMask)
        model.getImage().getArray()[:,:] = array
        residuals = lsst.afw.image.MaskedImageF(data, True)
        residuals -= model
        mosaic = lsst.afw.display.utils.Mosaic()
        mosaic.setMode("x")
        mosaic.append(data, "data")
        mosaic.append(model, "model")
        mosaic.append(residuals, "data-model")
        grid = mosaic.makeMosaic()
        lsst.afw.display.ds9.mtv(grid)
        lsst.afw.display.ds9.setMaskTransparency(85)
        for i in range(3):
            ds9box = mosaic.getBBox(i)
            center = ellipse.getCenter() + lsst.afw.geom.Extent2D(ds9box.getMin() - bbox.getMin())
            lsst.afw.display.ds9.dot(ellipse.getCore(), center.getX(), center.getY())

    def displayExposure(self, frame=0, doLabel=False):
        lsst.afw.display.ds9.mtv(self.exposure, frame=frame)
        ellipseKey = self.modelfits.schema.find("ref.ellipse").key
        centerKey = self.modelfits.schema.find("ref.center").key
        with lsst.afw.display.ds9.Buffering():
            for record in self.modelfits:
                lsst.afw.display.ds9.dot(record.get(ellipseKey),
                                         record.get(centerKey.getX()), record.get(centerKey.getY()),
                                         ctype="green")
                if doLabel:
                    lsst.afw.display.ds9.dot(str(record.getId()),
                                             record.get(centerKey.getX()), record.get(centerKey.getY()),
                                             ctype="cyan")
