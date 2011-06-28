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
import lsst.afw.display.utils
import lsst.afw.display.ds9 as ds9
import lsst.afw.geom.ellipses
import lsst.meas.multifit
import numpy
from matplotlib import pyplot
from lsst.meas.multifitData import DatasetMapper
from lsst.daf.persistence import ButlerFactory
from lsst.pex.policy import Policy
import lsst.meas.algorithms

class Viewer(object):

    def __init__(self, dataset, basisSize=8, psfShapeletOrder=4, nTestPoints=5, nGrowFp=1, 
                usePixelWeights=False, fitDeltaFunction=True, 
                isEllipticityActive=True, isRadiusActive=True, isPositionActive=False,
                maskPlaneNames=["BAD", "INTRP", "EDGE", "CR", "SAT"]):
        self.basis = lsst.meas.multifit.SourceMeasurement.loadBasis(basisSize)
        
        self.sourceMeasurement = lsst.meas.multifit.SourceMeasurement(basisSize, psfShapeletOrder,
                nTestPoints, nGrowFp, usePixelWeights, fitDeltaFunction, 
                isEllipticityActive, isRadiusActive, isPositionActive,
                maskPlaneNames)

        bf = ButlerFactory(mapper=DatasetMapper())
        butler = bf.create()
        self.psf = butler.get("psf", id=dataset)

        self.exposure = butler.get("exp", id=dataset)
        self.exposure.setPsf(self.psf)
        self.sources = butler.get("src", id=dataset)
        self.bitmask = lsst.afw.image.MaskU.getPlaneBitMask(maskPlaneNames)
        self.scaleFactor = 5
        self.fits = {}
        self.plots = {}
    
        self.m = lsst.afw.display.utils.Mosaic()
        self.gutter = 10
        self.m.setGutter(self.gutter)
        self.m.setBackground(0)
        self.m.setMode("x")

    def fit(self, index):
        if index in self.fits:
            return

        #fit the source 
        source = self.sources[index]
        self.sourceMeasurement.measure(self.exposure, source)

        #get fit parameters/coeff
        fitFp = self.sourceMeasurement.getFootprint()
        ellipse = self.sourceMeasurement.getEllipse()
        core = ellipse.getCore()
        coeff = self.sourceMeasurement.getCoeff()
        param = self.sourceMeasurement.getParam()
        #generate a larger, unmasked fp
        bbox = fitFp.getBBox()
        bounds = lsst.afw.geom.ellipses.Ellipse(ellipse)
        bounds.scale(self.scaleFactor)
        bbox.include(lsst.afw.geom.Box2I(bounds.computeEnvelope()))
        bbox.grow(10)
        bbox.clip(self.exposure.getBBox(lsst.afw.image.PARENT))
        fp = lsst.afw.detection.Footprint(bbox)       



        #generate a new evalution of the fitted model over this larger fp
        definition = self.sourceMeasurement.getEvaluator().getGrid().makeDefinition()
        dataArray = numpy.zeros(fp.getArea(), dtype=numpy.float32)

        lsst.afw.detection.flattenArray(
            fp,
            self.exposure.getMaskedImage().getImage().getArray(), 
            dataArray, 
            self.exposure.getXY0()
        )

        definition.frames[0].setFootprint(fp)
        definition.frames[0].setData(dataArray)

        if self.sourceMeasurement.usePixelWeights():
            weightArray = numpy.ones_like(dataArray)
            lsst.afw.detection.flattenArray(
                fp,
                self.exposure.getMaskedImage().getVariance().getArray(), 
                weightArray, 
                self.exposure.getXY0()
            )
            definition.frames[0].setWeights(weightArray)
        
        evaluator = lsst.meas.multifit.Evaluator.make(definition)
        evaluation = lsst.meas.multifit.Evaluation(evaluator)
        evaluation.update(
            param,
            coeff
        )        


        #construct MaskedImages of the data, model, and residuals
        #of the fitted model evaluated over the larger fp                        

        original = self.exposure.getMaskedImage().getImage().Factory(
                self.exposure.getMaskedImage().getImage(), 
                bbox, lsst.afw.image.PARENT)        
        img = lsst.afw.image.ImageD(bbox)        
        msk = lsst.afw.image.MaskU(bbox)        
        img.getArray()[:,:] = original.getArray().astype(numpy.float)[:,:]
        ds9.mtv(img, frame = 5)

        ones = numpy.ones(fitFp.getArea(), dtype=int)
        lsst.afw.detection.expandArray(fitFp, ones, msk.getArray(), bbox.getMin())
        mi = lsst.afw.image.makeMaskedImage(img, msk)

        modelArray = evaluation.getModelVector().astype(img.getArray().dtype)
        residualArray = evaluation.getResiduals().astype(img.getArray().dtype)

        model = img.Factory(bbox)
        lsst.afw.detection.expandArray(fp, modelArray, model.getArray(), bbox.getMin())
        ds9.mtv(model, frame = 3)
        modelMi = mi.Factory(model, msk)


        residual = img.Factory(bbox)
        lsst.afw.detection.expandArray(fp, residualArray, residual.getArray(), bbox.getMin())    
        ds9.mtv(residual, frame= 4)
        residualMi = mi.Factory(residual, msk)

      
        psfImg = lsst.afw.image.ImageD(bbox)
        psfModel = lsst.afw.image.ImageD(bbox)
        psfResidual = lsst.afw.image.ImageD(bbox)

        psfDataArray = numpy.zeros(fp.getArea(), float)
        psfModelArray = numpy.zeros(fp.getArea(), float)

        source = evaluator.getGrid().sources[0]
        localPsfData = source.getLocalPsf()
        shapelet = localPsfData.computeShapelet(
            lsst.afw.math.shapelets.HERMITE,
            self.sourceMeasurement.getPsfShapeletOrder()
        )
        multiShapelet = lsst.afw.math.shapelets.MultiShapeletFunction(shapelet)
        localPsfModel = lsst.afw.detection.ShapeletLocalPsf(localPsfData.getPoint(), multiShapelet)
        localPsfData.evaluatePointSource(fp, psfDataArray)
        localPsfModel.evaluatePointSource(fp, psfModelArray)
        psfResidualArray = psfDataArray - psfModelArray

        lsst.afw.detection.expandArray(fp, psfDataArray, psfImg.getArray(), bbox.getMin())
        lsst.afw.detection.expandArray(fp, psfModelArray, psfModel.getArray(), bbox.getMin())
        lsst.afw.detection.expandArray(fp, psfResidualArray, psfResidual.getArray(), bbox.getMin())
    
        d = {}
        d["fp"] = fitFp
        d["bbox"] = bbox
        d["ellipse"] = ellipse
        d["evaluator"] = evaluator
        d["evaluation"] = evaluation
        d["images"] = [mi, modelMi, residualMi]
        d["image_labels"] = ["image", "model", "residuals"]
        d["psf_images"] = [psfImg, psfModel, psfResidual]
        d["psf_labels"] = ["psf", "ShapeletLocalPsf model", "residuals"]
        self.fits[index] = d

    def plot(self, index):
        self.fit(index)
        d = self.fits[index]
        bbox = d["bbox"]
        ellipse = d["ellipse"]
        q = lsst.afw.geom.ellipses.Quadrupole(ellipse.getCore())

        center = ellipse.getCenter() - lsst.afw.geom.ExtentD(bbox.getMin())

        mosaic = self.m.makeMosaic(d["images"], frame=0)
        self.m.drawLabels(d["image_labels"], frame=0)
        
        ds9.dot("x", c=center.getX(), r=center.getY(), frame=0)
        ds9.dot("@:%f,%f,%f"%(q.getIXX(), q.getIXY(), q.getIYY()), center.getX(), center.getY(), frame=0)
        center += lsst.afw.geom.ExtentD(bbox.getDimensions() + lsst.afw.geom.ExtentI(self.gutter - 2))
        ds9.dot("x", c=center.getX(), r=center.getY(), frame=0)
        ds9.dot("@:%f,%f,%f"%(q.getIXX(), q.getIXY(), q.getIYY()), center.getX(), center.getY(), frame=0)
        center += lsst.afw.geom.ExtentD(bbox.getDimensions() + lsst.afw.geom.ExtentI(self.gutter - 2))
        ds9.dot("x", c=center.getX(), r=center.getY(), frame=0)
        ds9.dot("@:%f,%f,%f"%(q.getIXX(), q.getIXY(), q.getIYY()), center.getX(), center.getY(), frame=0)
        
        psfMosaic = self.m.makeMosaic(d["psf_images"], frame=1)
        self.m.drawLabels(d["psf_labels"], frame=1)        

