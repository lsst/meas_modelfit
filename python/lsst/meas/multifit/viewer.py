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
    dataVector = evaluation.getEvaluator().getDataVector()
    modelVector = evaluation.getModelVector()
    residualVector = evaluation.getResiduals()
    for source in grid.sources:
        pyplot.figure()
        pyplot.title("Object %d, Frame %d" % (source.object.id, source.frame.id))
        pixelIndices = slice(source.frame.getPixelOffset(), source.frame.getPixelCount())
        dataSubset = dataVector[pixelIndices]
        modelSubset = modelVector[pixelIndices]
        residualSubset = residualVector[pixelIndices]
        footprint = source.frame.getFootprint()
        ibox = footprint.getBBox()
        dbox = lsst.afw.geom.BoxD(ibox)
        extent = (dbox.getMinX(), dbox.getMaxX(), dbox.getMinY(), dbox.getMaxY())
        images = numpy.zeros((3, ibox.getHeight(), ibox.getWidth()), dtype=float)
        lsst.afw.detection.expandArray(footprint, dataSubset, images[0], bbox.getMin());
        lsst.afw.detection.expandArray(footprint, modelSubset, images[1], bbox.getMin());
        lsst.afw.detection.expandArray(footprint, residualSubset, image[2], bbox.getMin());
        vmin = images.min()
        vmax = images.max()
        for i in range(3):
            pyplot.subplot(1, 3, i + 1)
            pyplot.imshow(images[i], origin='lower', interpolation='nearest', 
                          vmin=vmin, vmax=vmax, extent=extent)
            if source.object.getBasis():
                ellipse = source.object.makeEllipse(evaluation.getParameters())
                ellipse.plot(fill=False)
            else:
                point = source.object.makePoint(evaluation.getParameters())
                pyplot.plot([self.point.getX()], [self.point.getY()], 'kx')
    pyplot.show()

def plotInterpreter(interpreter):
    evaluator = lsst.meas.multifit.Evaluator(interpreter.getGrid())
    evaluation = lsst.meas.multifit.Evaluation(evaluator, intepreter.computeParameterMean());
    evaluation.setCoefficients(intepreter.computeCoefficientMean())
    plotEvaluation(evaluation, intepreter.getGrid())

class Viewer(object):

    def __init__(self, dataset):
        self.policy = Policy()
        self.policy.add("nGrowFp", 3)
        self.policy.add("isVariable", False)
        self.policy.add("isPositionActive", False)
        self.policy.add("isRadiusActive", True)
        self.policy.add("isEllipticityActive", True)
        self.policy.add("maskPlaneName", "BAD")
        self.policy.add("basisName", "ed+15:4000")
        bf = ButlerFactory(mapper=DatasetMapper())
        butler = bf.create()
        self.psf = butler.get("psf", id=dataset)
        print self.psf
        self.exposure = butler.get("exp", id=dataset)
        self.exposure.setPsf(self.psf)
        self.sources = butler.get("src", id=dataset)
        self.bitmask = utils.makeBitmask(self.exposure.getMaskedImage().getMask(),
                                         self.policy.getArray("maskPlaneName"))
    
    def plot(self, source):
        interpreter = utils.fitSource(self.exposure, source, self.bitmask, self.policy)
