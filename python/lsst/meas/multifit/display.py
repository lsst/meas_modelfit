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
import numpy
from matplotlib import pyplot
import mpl_toolkits.mplot3d

import lsst.pex.config
import lsst.daf.persistence
import lsst.meas.extensions.multiShapelet
import lsst.afw.image
import lsst.afw.geom.ellipses
import lsst.afw.display.ds9

from .measureCcd import *
from .measureCoadd import *
from .multifitLib import *

__all__ = ("makeHistogramGrid", "InteractiveFitter")

def makeHistogramGrid(figure, names):
    """Function to create a grid of subplots for use in plotting
    the an N-dimensional density.

    The returned grid is a object-type numpy array containing
    matplotlib.axes.Axes objects.  The off-diagonal subplots are
    intended for contour plots, and the diagonal plots are
    intended for 1-d histogram plots.
    """
    n = len(names)
    grid = numpy.empty((n,n), dtype=object)
    for j in range(n):
        grid[0,j] = figure.add_subplot(n, n, 1 + j)
        grid[0,j].set_title(names[j])
    for i in range(1, n):
        grid[i,0] = figure.add_subplot(n, n, 1 + i*n, sharex=grid[0,0])
        grid[i,0].set_ylabel(names[i])
    for i in range(1, n):
        for j in range(1, n):
            grid[i,j] = figure.add_subplot(n, n, 1 + i*n + j, sharex=grid[0,j], sharey=grid[i,0])
    for i in range(n-1):
        for j in range(n):
            for tl in grid[i,j].get_xticklabels(): tl.set_visible(False)
    for i in range(n):
        for j in range(1,n):
            for tl in grid[i,j].get_yticklabels(): tl.set_visible(False)
    for i in range(n):
        grid[i,i] = grid[i,i].twinx()
        for tl in grid[i,i].get_yticklabels():
            tl.set_color('r')
    return grid

class Interactive(object):
    """Interactive analysis helper class

    This class manages a butler, calexp, modelfits catalog, and an instance
    of the MeasureCcdTask, allowing individual objects to be re-fit and plotted.
    """

    def __init__(self, rerun, config=None, dataId=None, coadd=False):
        """Construct an interactive analysis object.

        @param[in]  rerun    Output directory, relative to $S13_DATA_DIR/output.
                             measureCcd.py (or measureCoadd.py if coadd=True) must
                             have been run (possibly with prepOnly=True) previously
                             with this output directory.
        @param[in]  config   MeasureCcdTask.ConfigClass instance; if None, it
                             will be loaded from disk.
        @param[in]  dataId   Butler data ID of the image to analyze.
        @param[in]  coadd    Whether to process on the coadd instead of single-frame.
        """
        if dataId is None:
            if not coadd:
                dataId = dict(visit=100, raft="2,2", sensor="1,1")
            else:
                dataId = dict(tract=0, patch="2,2", filter="i")
        root = os.environ["S13_DATA_DIR"]
        self.butler = lsst.daf.persistence.Butler(os.path.join(root, "output", rerun))
        if not coadd:
            self.dataRef = self.butler.dataRef("calexp", dataId=dataId)
            if config is None:
                config = self.butler.get("measureCcd_config", immediate=True)
            self.exposure = self.dataRef.get("calexp", immediate=True)
            self.modelfits = self.dataRef.get("modelfits", immediate=True)
            self.task = MeasureCcdTask(config=config)
        else:
            self.dataRef = self.butler.dataRef("deepCoadd_calexp", dataId=dataId)
            if config is None:
                config = self.butler.get("deep_measureCoadd_config", immediate=True)
            self.exposure = self.dataRef.get("deepCoadd_calexp", immediate=True)
            self.modelfits = self.dataRef.get("deepCoadd_modelfits", immediate=True)
            self.task = MeasureCoaddTask(config=config)

    def fit(self, index=0, id=None):
        """Re-fit the object indicated by the given record sequential index
        or source ID, returning the Struct object returned by
        MeasureImageTask.processObject.
        """
        if id is not None:
            record = self.modelfits.find(id)
        else:
            record = self.modelfits[index]
        return self.task.processObject(self.exposure, record)

    def plotSamples(self, record, threshold=1E-15, **kwds):
        """Plot the samples corresponding to the given ModelFitRecord
        as a 3-d scatter plot.

        The 'threshold' parameter is used to restrict which sample points
        are plotted; points with weight < max(weight)*threshold are skipped.

        Additional keyword arguments are passed to the matplotlib scatter
        command.
        """
        samples = record.getSamples()
        mean = samples.computeMean()
        median = samples.computeQuantiles(numpy.array([0.5]))
        if self.task.config.useRefCat:
            key = self.task.keys["ref.ellipse"]
            ellipse = lsst.afw.geom.ellipses.BaseCore.make(samples.getEllipseType(), record.get(key))
            truth = ellipse.getParameterVector()
        cat = samples.getCatalog().copy(deep=True)
        maxR = cat["parameters"][:,2].max()
        mask = cat["weight"] / cat["weight"].max() > threshold
        cat = cat[mask].copy()
        fig = pyplot.figure()
        ax = fig.add_subplot(111, projection='3d')
        kwds.setdefault('edgecolor', 'none')
        kwds.setdefault('cmap', pyplot.cm.Blues)
        kwds.setdefault('s', 16)
        patches = ax.scatter(cat["parameters"][:,0], cat["parameters"][:,1], cat["parameters"][:,2],
                             c=numpy.log10(cat["weight"]), **kwds)
        cax = pyplot.colorbar(patches)
        cax.set_label("log10(probability) + [arbitrary constant]")
        ax.plot(mean[0:1], mean[1:2], mean[2:3], "gx", label="mean", markersize=9, markeredgewidth=2)
        ax.plot(median[0:1], median[1:2], median[2:3], "r+", label="median", markersize=9, markeredgewidth=2)
        if self.task.config.useRefCat:
            ax.plot(truth[0:1], truth[1:2], truth[2:3], "m^", label="truth", markersize=9, markeredgewidth=0)
        ax.legend()
        ax.set_xlabel("e1")
        ax.set_xlim(-1, 1)
        ax.set_ylabel("e2")
        ax.set_ylim(-1, 1)
        ax.set_zlabel("radius")
        ax.set_zlim(0, maxR)
        pyplot.show()
        return ax

    def displayResiduals(self, r, parameters=None):
        """Display the data postage stamp along with the model image and residuals in ds9.

        @param[in] r            Result Struct returned by fit()
        @param[in] parameters   Parameter vector for the model; defaults to the posterior mean.
        """
        center = r.record.getPointD(self.task.keys["source.center"])
        samples = r.record.getSamples()
        if parameters == "truth":
            key = self.task.keys["ref.ellipse"]
            truth = lsst.afw.geom.ellipses.BaseCore.make(samples.getEllipseType(), r.record.get(key))
            ellipse = lsst.afw.geom.ellipses.Ellipse(truth, center)
        else:
            if parameters == "mean" or parameters is None:
                parameters = samples.computeMean()
            elif parameters == "median":
                parameters = samples.computeQuantiles(numpy.array([0.5]))
            else:
                parameters = numpy.array(parameters)
            ellipse = r.record.getSamples().interpret(parameters, center)
        joint = r.objective.evaluate(ellipse)
        bbox = r.record.getFootprint().getBBox()
        bbox.grow(2)
        xe = numpy.arange(bbox.getBeginX(), bbox.getEndX(), dtype=numpy.float32)
        ye = numpy.arange(bbox.getBeginY(), bbox.getEndY(), dtype=numpy.float32)
        xg, yg = numpy.meshgrid(xe, ye)
        matrixBuilder = lsst.shapelet.MultiShapeletMatrixBuilderF(
            self.task.basis, r.psf, xg.ravel(), yg.ravel(), self.task.config.objective.useApproximateExp
            )
        matrix = numpy.zeros((joint.mu.size, xg.size), dtype=numpy.float32).transpose()
        matrixBuilder.build(matrix, ellipse)
        array = numpy.dot(matrix, joint.mu).reshape(*xg.shape)

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
