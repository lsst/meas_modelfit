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

from ..measureCcd import MeasureCcdTask
from ..measureCoadd import MeasureCoaddTask
from ..measureMulti import MeasureMultiTask
from .. import multifitLib

from .densityPlot import *
from .modelFitAdapters import *

__all__ = ("Interactive", )

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
            self.task = MeasureCcdTask(config=config)
        elif self.mode == "coadd":
            if dataId is None:
                dataId = dict(tract=0, patch="2,2", filter="i")
            self.dataRef = self.butler.dataRef("deepCoadd_calexp", dataId=dataId)
            if config is None:
                config = self.butler.get("deep_measureCoadd_config", immediate=True)
            self.task = MeasureCoaddTask(config=config)
        elif self.mode.startswith("multi"):
            if dataId is None:
                dataId = dict(tract=0, patch="2,2", filter="i")
            self.dataRef = self.butler.dataRef("deepCoadd_calexp", dataId=dataId)
            if config is None:
                config = self.butler.get("deep_measureMulti_config", immediate=True)
            self.task = MeasureMultiTask(config=config)
        self.inputs = self.task.readInputs(self.dataRef)
        self.cat = self.task.prepCatalog(self.inputs)

    def fit(self, index=0, id=None):
        """Re-fit the object indicated by the given record sequential index or source ID,
        returning the record.
        """
        if id is not None:
            record = self.cat.find(id)
        else:
            record = self.cat[index]
        likelihood = self.task.makeLikelihood(self.inputs, record)
        self.task.fitter.run(likelihood, record)
        return record

    def plotSamples(self, record):
        """Plot the samples and proposal distribution from a ModelFitRecord.
        """
        data = ModelFitDataAdapter(record)
        figure = matplotlib.pyplot.figure(record.getId(), figsize=(10, 10))
        p = DensityPlot(figure, data)
        p.layers["samples"] = HistogramLayer()
        p.layers["proposal"] = SurfaceLayer()
        p.draw()
        return p

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
        joint = r.likelihood.evaluate(ellipse)
        bbox = r.record.getFootprint().getBBox()
        bbox.grow(2)
        xe = numpy.arange(bbox.getBeginX(), bbox.getEndX(), dtype=numpy.float32)
        ye = numpy.arange(bbox.getBeginY(), bbox.getEndY(), dtype=numpy.float32)
        xg, yg = numpy.meshgrid(xe, ye)
        matrixBuilder = lsst.shapelet.MultiShapeletMatrixBuilderF(
            self.task.basis, r.psf, xg.ravel(), yg.ravel(), self.task.config.likelihood.useApproximateExp
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
