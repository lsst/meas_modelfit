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
        self.mode = mode.lower()
        if self.mode.startswith("test"):
            self.dataRef = rerun
            if config is None:
                config = MeasureCcdTask.ConfigClass()
            self.task = MeasureCcdTask(config=config)
        else:
            root = os.environ["S13_DATA_DIR"]
            self.butler = lsst.daf.persistence.Butler(os.path.join(root, "output", rerun))
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

    def displayResiduals(self, record, nonlinear="ref", amplitudes="ref"):
        """Display the data postage stamp along with the model image and residuals in ds9.

        @param[in] record       ModelFitRecord defining the object to display
        @param[in] nonlinear    Vector of nonlinear parameters, or a string tag (see below)
        @param[in] amplitudes   Vector of linear parameters, or a string tag (see below)

        String tags for the parameter arguments may be:
          mean -- use record.getIntepreter().compute*Mean()
          ref --- use record.get("ref.*")
          lstsq - solve the (unweighted) linear least-squares problem and use the solution
        """
        likelihood = self.task.makeLikelihood(self.inputs, record)

        if nonlinear == "mean":
            nonlinear = record.getInterpreter().computeNonlinearMean(record)
        elif nonlinear == "ref":
            nonlinear = record.get("ref.nonlinear")
        else:
            assert nonlinear.shape == (likelihood.getNonlinearDim(),)

        matrix = numpy.zeros((likelihood.getAmplitudeDim(), likelihood.getDataDim()),
                             dtype=multifitLib.Pixel).transpose()
        likelihood.computeModelMatrix(matrix, nonlinear)

        if amplitudes == "mean":
            amplitudes = record.getInterpreter().computeAmplitudeMean()
        elif amplitudes == "ref":
            amplitudes = record.get("ref.amplitudes")
        elif amplitudes == "lstsq":
            amplitudes, chisq, rank, sv = numpy.linalg.lstsq(matrix, likelihood.getData())
        else:
            assert amplitudes.shape == (likelihood.getAmplitudeDim(),)

        bbox = record.getFootprint().getBBox()
        bbox.grow(2)
        flatModelW = numpy.zeros(likelihood.getDataDim(), dtype=multifitLib.Pixel)
        flatModelW[:] = numpy.dot(matrix, amplitudes)
        flatDataW = likelihood.getData()

        imgDataU = lsst.afw.image.MaskedImageF(self.inputs.exposure.getMaskedImage(), bbox,
                                               lsst.afw.image.PARENT, True)
        bitmask = imgDataU.getMask().addMaskPlane("FIT_REGION")
        regionMask = lsst.afw.image.MaskU(bbox)
        lsst.afw.detection.setMaskFromFootprint(regionMask, record.getFootprint(), bitmask)
        dataMask = imgDataU.getMask()
        dataMask |= regionMask
        imgDataW = imgDataU.clone()
        imgDataW.getImage().set(0.0)
        imgDataW.getVariance().set(0.0)
        lsst.afw.detection.expandArray(record.getFootprint(), flatDataW, imgDataW.getImage().getArray(),
                                       imgDataW.getXY0())
        imgModelW = lsst.afw.image.MaskedImageF(lsst.afw.image.ImageF(bbox), regionMask)
        lsst.afw.detection.expandArray(record.getFootprint(), flatModelW, imgModelW.getImage().getArray(),
                                       imgModelW.getXY0())
        imgResidualsW = lsst.afw.image.MaskedImageF(imgDataW, True)
        imgResidualsW -= imgModelW
        mosaic = lsst.afw.display.utils.Mosaic()
        mosaic.setMode("x")
        mosaic.append(imgDataW, "data")
        mosaic.append(imgModelW, "model")
        mosaic.append(imgResidualsW, "data-model")
        grid = mosaic.makeMosaic()
        lsst.afw.display.ds9.mtv(grid)
        lsst.afw.display.ds9.setMaskTransparency(85)
