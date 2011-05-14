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

def plotEvaluation(evaluation, grid):
    print evaluation.getCoefficients()
    dataVector = evaluation.getEvaluator().getDataVector()
    print "dataVector.sum():", dataVector.sum()
    modelMatrix = evaluation.getModelMatrix()
    modelVector = evaluation.getModelVector()
    print "modelVector.sum():", modelVector.sum()
    residualVector = evaluation.getResiduals()
    for source in grid.sources:
        pyplot.figure()
        pyplot.title("Object %d, Frame %d" % (source.object.id, source.frame.id))
        pixelIndices = slice(source.frame.getPixelOffset(), source.frame.getPixelCount())
        dataSubset = dataVector[pixelIndices]
        modelSubset = modelVector[pixelIndices]
        residualSubset = residualVector[pixelIndices]
        footprint = source.frame.getFootprint()
        psfDataVector = numpy.zeros(footprint.getArea(), dtype=float)
        psfModelVector = numpy.zeros(footprint.getArea(), dtype=float)
        localPsfData = source.getLocalPsf()
        shapelet = localPsfData.computeShapelet(lsst.afw.math.shapelets.HERMITE,
                                                lsst.meas.multifit.ShapeletModelBasis.getPsfShapeletOrder())
        multiShapelet = lsst.afw.math.shapelets.MultiShapeletFunction(shapelet)
        localPsfModel = lsst.afw.detection.ShapeletLocalPsf(localPsfData.getPoint(), multiShapelet)
        localPsfData.evaluatePointSource(footprint, psfDataVector)
        localPsfModel.evaluatePointSource(footprint, psfModelVector)
        psfResidualVector = psfModelVector - psfDataVector
        ibox = footprint.getBBox()
        dbox = lsst.afw.geom.BoxD(ibox)
        extent = (dbox.getMinX(), dbox.getMaxX(), dbox.getMinY(), dbox.getMaxY())
        images = numpy.zeros((6, ibox.getHeight(), ibox.getWidth()), dtype=float)
        lsst.afw.detection.expandArray(footprint, dataSubset, images[0], ibox.getMin());
        lsst.afw.detection.expandArray(footprint, modelSubset, images[1], ibox.getMin());
        lsst.afw.detection.expandArray(footprint, residualSubset, images[2], ibox.getMin());
        lsst.afw.detection.expandArray(footprint, psfDataVector, images[3], ibox.getMin());
        lsst.afw.detection.expandArray(footprint, psfModelVector, images[4], ibox.getMin());
        lsst.afw.detection.expandArray(footprint, psfResidualVector, images[5], ibox.getMin());
        vmin = images[0].min()
        vmax = images[0].max()
        for i in range(3):
            pyplot.subplot(2, 3, i + 1)
            pyplot.imshow(images[i], origin='lower', interpolation='nearest', 
                          vmin=vmin, vmax=vmax, extent=extent)
            if source.object.getBasis():
                ellipse = source.object.makeEllipse(evaluation.getParameters())
                ellipse.plot(fill=False)
            else:
                point = source.object.makePoint(evaluation.getParameters())
                pyplot.plot([self.point.getX()], [self.point.getY()], 'kx')
            print i, images[i].sum()
        vmin = images[3].min()
        vmax = images[3].max()
        psfEllipse = localPsfData.computeMoments()
        for i in range(3):
            pyplot.subplot(2, 3, i + 4)
            pyplot.imshow(images[i + 3], origin='lower', interpolation='nearest', 
                          vmin=vmin, vmax=vmax, extent=extent)
            psfEllipse.plot(fill=False)
            print i+3, images[i + 3].sum()
    pyplot.show()

def plotInterpreter(interpreter):
    evaluator = lsst.meas.multifit.Evaluator.make(interpreter.getGrid())
    if evaluator.getParameterSize() > 0:
        evaluation = lsst.meas.multifit.Evaluation(evaluator, interpreter.computeParameterMean())
    else:
        evaluation = lsst.meas.multifit.Evaluation(evaluator)
    coefficients = interpreter.computeCoefficientMean()
    evaluation.setCoefficients(coefficients)
    plotEvaluation(evaluation, interpreter.getGrid())

def makeBitMask(mask, maskPlaneNames):
    bitmask=0
    for name in maskPlaneNames:
        bitmask |= mask.getPlaneBitMask(name)    
    return bitmask

def makeGrid(exposure, source, policy):
    bitmask = makeBitMask(exposure.getMaskedImage().getMask(), policy.getArray("maskPlaneName"))
    quadrupole = lsst.afw.geom.ellipses.Quadrupole(source.getIxx(), source.getIyy(), source.getIxy())
    ellipse = lsst.afw.geom.ellipses.Ellipse(
        lsst.meas.multifit.EllipseCore(quadrupole), 
        lsst.afw.geom.Point2D(source.getXAstrom(), source.getYAstrom())
        )
    basis = utils.loadBasis(policy.get("basisName"))
    footprint = lsst.afw.detection.growFootprint(source.getFootprint(), policy.get("nGrowFp"))
    definition = lsst.meas.multifit.Definition.make(
        exposure, footprint, basis, ellipse,
        policy.get("isEllipticityActive"),
        policy.get("isRadiusActive"),
        policy.get("isPositionActive"),
        bitmask
        )
    return lsst.meas.multifit.Grid.make(definition)

class Viewer(object):

    def __init__(self, dataset):
        self.policy = Policy()
        self.policy.add("nGrowFp", 3)
        self.policy.add("isVariable", False)
        self.policy.add("isPositionActive", False)
        self.policy.add("isRadiusActive", False)
        self.policy.add("isEllipticityActive", False)
        self.policy.add("maskPlaneName", "BAD")
        self.policy.add("basisName", "ed+06:2000")
        self.policy.add("ftol", 1E-3)
        self.policy.add("gtol", 1E-3)
        self.policy.add("minStep", 1E-8)
        self.policy.add("maxIter", 200)
        self.policy.add("tau", 1E-3)
        self.policy.add("retryWithSvd", True)
        bf = ButlerFactory(mapper=DatasetMapper())
        butler = bf.create()
        self.psf = butler.get("psf", id=dataset)
        self.exposure = butler.get("exp", id=dataset)
        self.exposure.setPsf(self.psf)
        self.sources = butler.get("src", id=dataset)
        self.bitmask = utils.makeBitMask(self.exposure.getMaskedImage().getMask(),
                                         self.policy.getArray("maskPlaneName"))
        self.optimizer = lsst.meas.multifit.GaussNewtonOptimizer()
    
    def plot(self, source):
        self.grid = makeGrid(self.exposure, source, self.policy)
        self.evaluator = lsst.meas.multifit.Evaluator.make(self.grid)
        self.distribution = self.optimizer.solve(
            self.evaluator,
            self.policy.get("ftol"),
            self.policy.get("gtol"),
            self.policy.get("minStep"),
            self.policy.get("maxIter"),
            self.policy.get("tau"),
            self.policy.get("retryWithSvd")
            )
        if self.distribution is None:
            print "Optimizer returned None"
            return
        interpreter = lsst.meas.multifit.UnifiedSimpleInterpreter.make(
            self.distribution,
            self.grid
            )
        print "model flux:", interpreter.computeFluxMean(0, 0)
        plotInterpreter(interpreter)
        #self.evaluation = lsst.meas.multifit.Evaluation(self.evaluator)
        #plotEvaluation(self.evaluation, self.grid)
