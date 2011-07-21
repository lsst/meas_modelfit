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

def makeViewFootprint(sourceMeasurement, exposure, bboxGrow=10, ellipseScale=5):
    fitFp = sourceMeasurement.getFootprint()
    ellipse = lsst.afw.geom.ellipses.Ellipse(sourceMeasurement.getEllipse())
    ellipse.scale(ellipseScale)
    bbox = fitFp.getBBox()
    bbox.include(lsst.afw.geom.BoxI(ellipse.computeEnvelope()))
    bbox.grow(10)
    bbox.clip(exposure.getBBox(lsst.afw.image.PARENT))

    return lsst.afw.detection.Footprint(bbox)


class Viewer(object):

    def __init__(self, dataset, **kw):
        self.sourceMeasurement, policy = lsst.meas.multifit.makeSourceMeasurement(**kw)

        bf = ButlerFactory(mapper=DatasetMapper())
        butler = bf.create()
        self.psf = butler.get("psf", id=dataset)

        self.exposure = butler.get("exp", id=dataset)
        self.exposure.setPsf(self.psf)
        self.sources = butler.get("src", id=dataset)
        self.bitmask = lsst.afw.image.MaskU.getPlaneBitMask(
            self.sourceMeasurement.getOptions().maskPlaneNames
            )
        self.scaleFactor = 5
        self.fits = {}
        self.plots = {}
    
        self.m = lsst.afw.display.utils.Mosaic()
        self.gutter = 3
        self.m.setGutter(self.gutter)
        self.m.setBackground(0)
        self.m.setMode("x")

    def fit(self, index):
        if index in self.fits:
            return

        #fit the source 
        source = self.sources[index]
        self.sourceMeasurement.measure(self.exposure, source)
        coefficients = self.sourceMeasurement.getCoefficients()
        parameters = self.sourceMeasurement.getParameters()
        fitFp = self.sourceMeasurement.getFootprint()

        #get fit parameters/coefficientks
        viewFp = makeViewFootprint(self.sourceMeasurement, self.exposure, 10, self.scaleFactor)
        bbox = viewFp.getBBox()

        #construct MaskedImages of the data, model, and residuals
        #of the fitted model evaluated over the larger fp                        
        sub = self.exposure.Factory(
                self.exposure, 
                bbox, lsst.afw.image.PARENT)        
        img = lsst.afw.image.ImageD(bbox)       
        msk = lsst.afw.image.MaskU(bbox)    
        msk.addMaskPlane("multifit")
        bitmask = msk.getPlaneBitMask("multifit")
        lsst.afw.detection.setMaskFromFootprint(msk, fitFp, bitmask)

        #get the inverted mask
        msk.getArray()[:,:] = numpy.logical_not(msk.getArray() == bitmask)*bitmask
        
        imgDtype = img.getArray().dtype
        img.getArray()[:,:] = sub.getMaskedImage().getImage().getArray().astype(imgDtype)[:,:]
        mi = lsst.afw.image.makeMaskedImage(img, msk)

        model = self.sourceMeasurement.getModelImage(sub)
        modelMi = mi.Factory(model, msk)

        residual = img.Factory(img, True)
        residual -= model
        residualMi = mi.Factory(residual, msk)
       
        #generate images of the psf and psf models
        psfImg = self.sourceMeasurement.getPsfImage(viewFp)
        psfModel = self.sourceMeasurement.getPsfModel(viewFp)
        psfResidual = psfImg.Factory(psfImg, True)
        psfResidual -= psfModel

        d = {}
        d["fp"] = fitFp
        d["bbox"] = bbox
        d["ellipse"] = self.sourceMeasurement.getEllipse()
        d["coefficients"] = coefficients
        d["evaluator"] = self.sourceMeasurement.getEvaluator()
        d["evaluation"] = self.sourceMeasurement.getEvaluation
        d["test_points"] = self.sourceMeasurement.getTestPoints()
        d["objective_value"] = self.sourceMeasurement.getObjectiveValue()
        d["images"] = [mi, modelMi, residualMi]
        d["image_labels"] = ["image", "model", "residuals"]
        d["psf_images"] = [psfImg, psfModel, psfResidual]
        d["psf_labels"] = ["psf", "ShapeletLocalPsf model", "residuals"]
        self.fits[index] = d

    def plot(self, index, mode="ds9"):
        self.fit(index)
        d = self.fits[index]
        bbox = d["bbox"]
        ellipse = d["ellipse"]
        q = lsst.afw.geom.ellipses.Quadrupole(ellipse.getCore())
        
        center = lsst.afw.geom.Point2D(ellipse.getCenter().getX() - bbox.getMinX(), ellipse.getCenter().getY() - bbox.getMinY())
        src = self.sources[index]

        if(mode=="ds9"):
            ds9.setMaskTransparency(75)

            #one frame for psf mosaic
            psfMosaic = self.m.makeMosaic(d["psf_images"], frame=1)
            self.m.drawLabels(d["psf_labels"], frame=1)  

            #another for the model mosaic
            mosaic = self.m.makeMosaic(d["images"], frame=0)
            self.m.drawLabels(d["image_labels"], frame=0)

            for frame in range(2):
                center = lsst.afw.geom.Point2D(
                    ellipse.getCenter().getX() - bbox.getMinX(), 
                    ellipse.getCenter().getY() - bbox.getMinY()
                )
                #plot initial, and fit ellipses over model mosaic
                #image
                #ds9.dot("x", c=center.getX(), r=center.getY(), frame=0)
                ds9.dot("@:%f,%f,%f"%(q.getIXX(), q.getIXY(), q.getIYY()), center.getX(), center.getY(), frame=frame)
                ds9.dot("@:%f,%f,%f"%(src.getIxx(), src.getIxy(), src.getIyy()), center.getX(), center.getY(), frame=frame, ctype=ds9.BLUE)
                #model
                center += lsst.afw.geom.ExtentD(bbox.getWidth() + self.gutter , 0)
                #ds9.dot("x", c=center.getX(), r=center.getY(), frame=0)
                ds9.dot("@:%f,%f,%f"%(q.getIXX(), q.getIXY(), q.getIYY()), center.getX(), center.getY(), frame=frame)
                ds9.dot("@:%f,%f,%f"%(src.getIxx(), src.getIxy(), src.getIyy()), center.getX(), center.getY(), frame=frame, ctype=ds9.BLUE)
                #residual
                center += lsst.afw.geom.ExtentD(bbox.getWidth() + self.gutter, 0)
                #ds9.dot("x", c=center.getX(), r=center.getY(), frame=0)
                ds9.dot("@:%f,%f,%f"%(q.getIXX(), q.getIXY(), q.getIYY()), center.getX(), center.getY(), frame=frame)
                ds9.dot("@:%f,%f,%f"%(src.getIxx(), src.getIxy(), src.getIyy()), center.getX(), center.getY(), frame=frame, ctype=ds9.BLUE)
            

            ds9.ds9Cmd("tile mode row")
            ds9.ds9Cmd("tile yes")

    def plotProfile(self, index):
        self.fit(index)
        d = self.fits[index]
        radii = numpy.linspace(0.0, 5.0, 200)
        fullProfile = numpy.zeros((radii.size, d['evaluator'].getCoefficientCount()), dtype=float)
        allProfiles = {"full":fullProfile, "psf": fullProfile.copy()}
        v = d["coefficients"] * self.sourceMeasurement.getIntegration()
        v /= v.sum()
        allFractions = {}
        offset = 0
        order = []
        if self.sourceMeasurement.getOptions().fitDeltaFunction:
            allFractions["psf"] = v[offset]
            order.append("psf")
            offset += 1
        if self.sourceMeasurement.getOptions().fitExponential:
            expProfile = numpy.zeros(fullProfile.shape, dtype=float)
            self.sourceMeasurement.getExponentialBasis().evaluateRadialProfile(
                expProfile[:,offset:offset+1], radii
                )
            allProfiles["exp"] = expProfile
            allFractions["exp"] = v[offset]
            order.append("exp")
            fullProfile += expProfile
            offset += 1
        if self.sourceMeasurement.getOptions().fitDeVaucouleur:
            devProfile = numpy.zeros(fullProfile.shape, dtype=float)
            self.sourceMeasurement.getDeVaucouleurBasis().evaluateRadialProfile(
                devProfile[:,offset:offset+1], radii
                )
            allProfiles["dev"] = devProfile
            allFractions["dev"] = v[offset]
            order.append("dev")
            fullProfile += devProfile
            offset += 1
        if self.sourceMeasurement.getOptions().shapeletOrder >= 0:
            shpProfile = numpy.zeros(fullProfile.shape, dtype=float)
            self.sourceMeasurement.getShapeletBasis().evaluateRadialProfile(
                shpProfile[:,offset:], radii
                )
            allProfiles["shp"] = shpProfile
            allFractions["shp"] = v[offset:].sum()
            order.append("shp")
            fullProfile += shpProfile
        f = d["ellipse"].getCore().getArea() / numpy.pi
        r = d["ellipse"].getCore().getTraceRadius()
        colors = {"full": "k", "exp":"b", "dev":"r", "shp":"g", "psf":"y"}
        pyplot.figure()
        ax1 = pyplot.axes([0.1, 0.38, 0.6, 0.24])
        ax1.plot(radii * r, numpy.dot(fullProfile, d["coefficients"]) / f, "k", label="combined")
        for k in order:
            if k not in allProfiles: continue
            v = allProfiles[k]
            ax1.plot(radii * r, numpy.dot(v, d["coefficients"]) / f, colors[k], label=k)
            pyplot.xlabel("radius (pixels)")
	pyplot.axvline(r, color="k", linestyle=":")
        pyplot.legend()
        ax2 = pyplot.axes([0.1, 0.66, 0.6, 0.24], sharex=ax1)
        ax2.semilogy(radii * r, numpy.dot(fullProfile, d["coefficients"]) / f, "k", label="combined")
        for k in order:
            if k not in allProfiles: continue
            v = allProfiles[k]
            ax2.semilogy(radii * r, numpy.dot(v, d["coefficients"]) / f, colors[k], label=k)
	pyplot.axvline(r, color="k", linestyle=":")
        pyplot.title("deconvolved radial profile")
        order.reverse()

        #objective value vs. radius
        ax3 = pyplot.axes([0.1, 0.1, 0.6, 0.24]) #, sharex=ax1)
        ax3.plot(d["test_points"][:, 2], d["objective_value"][:, 0, 0], "ok")
        ax3.axvline(r, color="k", linestyle=":")

        #flux fraction
        ax4 = pyplot.axes([0.75, 0.1, 0.2, 0.8])
        ax4.barh(numpy.arange(0, len(order)), [allFractions[k] for k in order], 
                 color=[colors[k] for k in order])
        pyplot.ylim(-0.2, len(order))
        pyplot.title("flux fraction")
        ax4.set_yticklabels([])
        pyplot.xlim(0, 1.0)

        

        pyplot.show()
