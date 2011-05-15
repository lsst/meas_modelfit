#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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

import lsst.afw.image
import lsst.afw.detection
import lsst.afw.image
import lsst.afw.geom.ellipses
import lsst.meas.multifit
import numpy
from . import utils
from matplotlib import pyplot
from lsst.meas.multifitData import DatasetMapper
from lsst.daf.persistence import ButlerFactory
from lsst.pex.policy import Policy
import lsst.meas.algorithms

def makeBitMask(mask, maskPlaneNames):
    bitmask=0
    for name in maskPlaneNames:
        bitmask |= mask.getPlaneBitMask(name)    
    return bitmask

def makeGrid(exposure, source, basis, bitmask, policy):
    quadrupole = lsst.afw.geom.ellipses.Quadrupole(source.getIxx(), source.getIyy(), source.getIxy())
    ellipse = lsst.afw.geom.ellipses.Ellipse(
        lsst.meas.multifit.EllipseCore(quadrupole), 
        lsst.afw.geom.Point2D(source.getXAstrom(), source.getYAstrom())
        )
    footprint = lsst.afw.detection.growFootprint(source.getFootprint(), policy.get("SHAPELET_MODEL_8.nGrowFp"))
    definition = lsst.meas.multifit.Definition.make(
        exposure, footprint, basis, ellipse,
        policy.get("SHAPELET_MODEL_8.isEllipticityActive"),
        policy.get("SHAPELET_MODEL_8.isRadiusActive"),
        policy.get("SHAPELET_MODEL_8.isPositionActive"),
        bitmask
        )
    return lsst.meas.multifit.Grid.make(definition)

class Viewer(object):

    def __init__(self, dataset):
        self.policy = Policy()
        self.policy.add("SHAPELET_MODEL_8.enabled", True)
        self.policy.add("SHAPELET_MODEL_8.nGrowFp", 3)
        self.policy.add("SHAPELET_MODEL_8.isPositionActive", False)
        self.policy.add("SHAPELET_MODEL_8.isRadiusActive", False)
        self.policy.add("SHAPELET_MODEL_8.isEllipticityActive", False)
        self.policy.add("SHAPELET_MODEL_8.maskPlaneName", "BAD")
        self.basis = utils.loadBasis("ed+06:2000")
        self.nTestPoints = 5

        bf = ButlerFactory(mapper=DatasetMapper())
        butler = bf.create()
        self.psf = butler.get("psf", id=dataset)

        expD = butler.get("exp", id=dataset)
        miD = expD.getMaskedImage()
        miF = lsst.afw.image.MaskedImageF(miD.getImage().convertF(), miD.getMask(), miD.getVariance())

        self.exposure = lsst.afw.image.ExposureF(miF, expD.getWcs())
        self.exposure.setPsf(self.psf)
        self.sources = butler.get("src", id=dataset)
        self.bitmask = utils.makeBitMask(self.exposure.getMaskedImage().getMask(),
                                         self.policy.getArray("SHAPELET_MODEL_8.maskPlaneName"))
        self.optimizer = lsst.meas.multifit.BruteForceSourceOptimizer()
        self.measurePhotometry = lsst.meas.algorithms.makeMeasurePhotometry(self.exposure)
        self.measurePhotometry.addAlgorithm("SHAPELET_MODEL_8")
        self.measurePhotometry.configure(self.policy)

        self.scaleFactor = 5
        self.fits = {}
        self.plots = {}
    
    def fit(self, index, mode=None):
        source = self.sources[index]
        d = {"grid": makeGrid(self.exposure, source, self.basis, self.bitmask, self.policy)}
        d["evaluator"] = lsst.meas.multifit.Evaluator.make(d["grid"])
        d["evaluation"] = lsst.meas.multifit.Evaluation(d["evaluator"])
        if mode == "usePhotometry":

            photom = self.measurePhotometry.measure(lsst.afw.detection.Peak(), source).find("SHAPELET_MODEL_8")

            nCoeff = d["evaluator"].getCoefficientSize()
            coefficients = numpy.zeros(nCoeff, dtype=float)
            for i in range(nCoeff):
                coefficients[i] = photom.get(i, lsst.afw.detection.Schema("COEFFICIENTS", 6, lsst.afw.detection.Schema.DOUBLE))
            d["evaluation"].setCoefficients(coefficients)
            d["flux"] = photom.getFlux()

        elif mode == "fitPython":
            self.optimizer.solve(d["evaluator"], self.nTestPoints)
            d["evaluation"].update(self.optimizer.getBestParameters(), self.optimizer.getCoefficients())
        self.fits[index] = d

    def makePlotDict(self, index):
        if index not in self.fits:
            self.fit(index)
        fitDict = self.fits[index]
        definition = fitDict["grid"].makeDefinition()
        for frame in definition.frames:
            bbox = frame.getFootprint().getBBox()
            for obj in fitDict["grid"].objects:
                bounds = obj.makeEllipse(fitDict["evaluation"].getParameters())
                bounds.scale(self.scaleFactor)
                bbox.include(lsst.afw.geom.Box2I(bounds.computeEnvelope()))
            bbox.clip(self.exposure.getBBox(lsst.afw.image.PARENT))
            frame.setFootprint(lsst.afw.detection.Footprint(bbox))
            dataArray = numpy.zeros(frame.getFootprint().getArea(), dtype=float)
            lsst.afw.detection.flattenArray(
                frame.getFootprint(),
                self.exposure.getMaskedImage().getImage().getArray(), 
                dataArray, 
                self.exposure.getXY0()
                )
            frame.setData(dataArray)
            frame.setWeights(numpy.ones(dataArray.shape, dtype=float))
        plotDict = {"grid": lsst.meas.multifit.Grid.make(definition)}
        plotDict["evaluator"] = lsst.meas.multifit.Evaluator.make(plotDict["grid"])
        plotDict["evaluation"] = lsst.meas.multifit.Evaluation(plotDict["evaluator"])
        plotDict["evaluation"].update(
            fitDict["evaluation"].getParameters(),
            fitDict["evaluation"].getCoefficients(),
            )
        self.plots[index] = plotDict

    def makeImages(self, d):
        dataVector = d["evaluation"].getEvaluator().getDataVector()
        modelVector = d["evaluation"].getModelVector()
        residualVector = d["evaluation"].getResiduals()
        d["frames"] = []
        ellipses = []
        points = []
        for obj in d["grid"].objects:
            if obj.getBasis():
                ellipses.append(obj.makeEllipse(d["evaluation"].getParameters()))
            else:
                points.append(obj.makePoint(d["evaluation"].getParameters()))
        for frame in d["grid"].frames:    
            pixelIndices = slice(frame.getPixelOffset(), frame.getPixelCount())
            dataSubset = dataVector[pixelIndices]
            modelSubset = modelVector[pixelIndices]
            residualSubset = residualVector[pixelIndices]
            footprint = frame.getFootprint()
            ibox = footprint.getBBox()
            dbox = lsst.afw.geom.BoxD(ibox)
            f = {}
            f["images"] = numpy.zeros((3, ibox.getHeight(), ibox.getWidth()), dtype=float)
            f["extent"] = (dbox.getMinX(), dbox.getMaxX(), dbox.getMinY(), dbox.getMaxY())
            lsst.afw.detection.expandArray(footprint, dataSubset, f["images"][0], ibox.getMin())
            lsst.afw.detection.expandArray(footprint, modelSubset, f["images"][1], ibox.getMin())
            lsst.afw.detection.expandArray(footprint, residualSubset, f["images"][2], ibox.getMin())
            f["vmin"] = f["images"][0].min()
            f["vmax"] = f["images"][0].max()
            f["ellipses"] = ellipses
            f["points"] = points
            d["frames"].append(f)

    def makePsfImages(self, d):
        d["psfs"] = []
        for obj in d["grid"].objects:
            l = []
            for source in obj.sources:
                footprint = source.frame.getFootprint()
                ibox = footprint.getBBox()
                dbox = lsst.afw.geom.BoxD(ibox)
                psfDataVector = numpy.zeros(footprint.getArea(), dtype=float)
                psfModelVector = numpy.zeros(footprint.getArea(), dtype=float)
                localPsfData = source.getLocalPsf()
                shapelet = localPsfData.computeShapelet(
                    lsst.afw.math.shapelets.HERMITE,
                    lsst.meas.multifit.ShapeletModelBasis.getPsfShapeletOrder()
                    )
                multiShapelet = lsst.afw.math.shapelets.MultiShapeletFunction(shapelet)
                localPsfModel = lsst.afw.detection.ShapeletLocalPsf(localPsfData.getPoint(), multiShapelet)
                localPsfData.evaluatePointSource(footprint, psfDataVector)
                localPsfModel.evaluatePointSource(footprint, psfModelVector)
                psfResidualVector = psfModelVector - psfDataVector
                f = {}
                f["images"] = numpy.zeros((3, ibox.getHeight(), ibox.getWidth()), dtype=float)
                f["extent"] = (dbox.getMinX(), dbox.getMaxX(), dbox.getMinY(), dbox.getMaxY())
                lsst.afw.detection.expandArray(footprint, psfDataVector, f['images'][0], ibox.getMin())
                lsst.afw.detection.expandArray(footprint, psfModelVector, f['images'][1], ibox.getMin())
                lsst.afw.detection.expandArray(footprint, psfResidualVector, f['images'][2], ibox.getMin())
                f["vmin"] = f["images"][0].min()
                f["vmax"] = f["images"][0].max()
                f["ellipses"] = [localPsfData.computeMoments()]
                f["points"] = []
                l.append(f)
            d["psfs"].append(l)

    def plotRow(self, frame, row, nrows):
        colors = ("r", "g", "b")
        for i in range(3):
            axes = pyplot.subplot(nrows, 4, 4 * row + 1 + i)
            pyplot.imshow(frame['images'][i], origin='lower', interpolation='nearest', 
                          extent=frame['extent'], vmin=frame['vmin'], vmax=frame['vmax'])
            axes.get_xaxis().set_ticks([])
            axes.get_yaxis().set_ticks([])
            for ellipse in frame['ellipses']:
                ellipse.plot(fill=False)
            for point in frame['points']:
                pyplot.plot([point.getX()], [point.getY()], 'kx')
            pyplot.xlabel("sum=%f" % frame["images"][i].sum())
            if i == 0:
                pyplot.ylabel("min=%f\nmax=%f" % (frame['vmin'], frame['vmax']))
            
            #hist, edges = numpy.histogram(frame["images"][i], bins=10, new=True)
            pyplot.subplot(nrows, 4, 4 * row + 4)
            #pyplot.plot(0.5 * (edges[:-1] + edges[1:]), hist, colors[i])
            pyplot.hist(frame["images"][i].ravel(), bins=20, alpha=0.2,
                        edgecolor=colors[i], facecolor=colors[i])

    def plot(self, index, mode="fitPython"):
        self.fit(index, mode=mode)
        self.makePlotDict(index)
        self.makeImages(self.fits[index])
        self.makeImages(self.plots[index])
        self.makePsfImages(self.fits[index])
        #print "Computed flux: ", self.fits[index]["flux"]
        pyplot.figure()
        for n in range(1):
            self.plotRow(self.fits[index]["frames"][n], 0, 3)
            self.plotRow(self.plots[index]["frames"][n], 1, 3)
            self.plotRow(self.fits[index]["psfs"][0][n], 2, 3)
        pyplot.show()
