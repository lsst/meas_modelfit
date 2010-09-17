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

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDet
import lsst.afw.math as afwMath
import lsst.afw.coord as afwCoord
import lsst.meas.multifit as measMult
import lsst.pex.policy as pexPolicy
import numpy
import numpy.random
import lsst.afw.display.ds9 as ds9
import lsst.pex.logging as pexLog


def makeModelExposure(model, psf, wcs, noiseFactor=0):
    fp = model.computeProjectionFootprint(psf, wcs)
    nPix = fp.getNpix()
    box = fp.getBBox()
    proj = model.makeProjection(psf, wcs, fp)
    modelImage = afwImage.MaskedImageF(box.getWidth(), box.getHeight())
    modelImage.setXY0(box.getX0(), box.getY0())
    
    imageVector = proj.computeModelImage()
    varianceVector = numpy.zeros(nPix, dtype=numpy.float32)

    measMult.expandImageF(fp, subImage, imageVector, varianceVector)

    afwRandom = afwMath.Random()
    randomImg = afwImage.ImageF(modelImage.getDimensions())
    afwMath.randomGaussianImage(randomImg, afwRandom)
    randomImg *= noiseFactor
    img = modelImage.getImage()
    img += randomImg

    stats = afwMath.makeStatistics(modelImage, afwMath.VARIANCE)
    variance = stats.getValue(afwMath.VARIANCE)
    modelImage.getVariance().set(variance)

    exp = afwImage.ExposureF(modelImage, wcs)
    exp.setPsf(psf)
    return exp

def applyFitter():
    #exp = afwImage.ExposureF("c00.fits")
    #wcs = exp.getWcs()
    frameId = 0
    crVal = afwGeom.makePointD(45,45)
    crPix = afwGeom.makePointD(1,1)
    wcs = afwImage.createWcs(crVal, crPix, 0.0001, 0., 0., 0.0001)
    print wcs.getCDMatrix()

    axes = afwGeom.ellipses.Axes(25,30,0)
    affine = wcs.linearizePixelToSky(crVal)
    print "affine", [affine[i] for i in range(6)]
    transformedAxes = axes.transform(affine)

    logShear = afwGeom.ellipses.LogShear(transformedAxes)
    print "logShear", logShear
    sersicIndex = 1.5
    flux = 35.0

    model = measMult.createSersicModel(flux, crVal, logShear, sersicIndex)

    psf = afwDet.createPsf("DoubleGaussian", 7, 7, 1.0)

    exp = makeModelExposure(model, psf, wcs, 20.0)
    
    ds9.mtv(exp, frame=frameId)
    frameId +=1

    expList = measMult.ExposureListF()
    expList.append(exp)
    
    jiggeredLogShear = afwGeom.ellipses.LogShear(logShear[0]*1.1, logShear[1]*1.1, logShear[2]*1.1)
    #flux *= 1.1
    testModel = measMult.createSersicModel(flux, crVal, jiggeredLogShear, sersicIndex)
    #testModel = model
    modelEvaluator = measMult.ModelEvaluator(testModel)
    modelEvaluator.setExposureList(expList)

    fitterPolicy = pexPolicy.Policy()

    fitter = measMult.MinuitNumericFitter(fitterPolicy)

    errors = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    result = fitter.apply(modelEvaluator, errors)

    
    #print "nIterations: %d"%result.sdqaMetrics.get("nIterations")
    #print "chisq: %d"%result.chisq
    #print "dChisq: %d"%result.dChisq
    print modelEvaluator.getLinearParameters()
    print modelEvaluator.getNonlinearParameters()
  
    print "nPix", modelEvaluator.getNPixels()
    print "nIterations", result.nIterations
    print "chisq", result.chisq

    om = modelEvaluator.getModel().clone()
    outExp = makeModelExposure(om, psf, wcs)

    ds9.mtv(oe, frame=frameId)



    
if __name__ == "__main__":
    applyFitter()
