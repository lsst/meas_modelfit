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


def applyFitter():
    #exp = afwImage.ExposureF("c00.fits")
    #wcs = exp.getWcs()
    crVal = afwGeom.makePointD(45,45)
    crPix = afwGeom.makePointD(1,1)
    wcs = afwImage.createWcs(crVal, crPix, 0.0001, 0., 0., 0.0001)
    print wcs.getCDMatrix()

    axes = afwGeom.ellipses.Axes(200,180,0)
    affine = wcs.linearizePixelToSky(crVal)
    print "affine", [affine[i] for i in range(6)]
    transformedAxes = axes.transform(affine)

    logShear = afwGeom.ellipses.LogShear(transformedAxes)
    print "logShear", logShear
    sersicIndex = 1.5
    flux = 10000.0

    model = measMult.createExponentialModel(flux, crVal, logShear)

    psf = afwDet.createPsf("DoubleGaussian", 2, 2, 1.0)

    fp = model.computeProjectionFootprint(psf, wcs)
    nPix = fp.getNpix()
    bbox = fp.getBBox()
    grownBbox = afwImage.BBox(bbox.getLLC(), bbox.getURC())
    grownBbox.shift(-30,-30)
    grownBbox.setWidth(bbox.getWidth() + 60)
    grownBbox.setHeight(bbox.getHeight() + 60)

    proj = model.makeProjection(psf, wcs, fp)
    modelImage = afwImage.MaskedImageF(grownBbox.getWidth(), grownBbox.getHeight())
    modelImage.setXY0(grownBbox.getX0(), grownBbox.getY0())

    imageVector = proj.computeModelImage()
    varianceVector = numpy.zeros(nPix, dtype=numpy.float32)

    bbox.shift(-bbox.getX0(), -bbox.getY0()) 
    subImage = modelImage.Factory(modelImage, bbox)
    measMult.expandImageF(fp, subImage, imageVector, varianceVector)
    
    afwRandom = afwMath.Random()
    randomImg = afwImage.ImageF(modelImage.getDimensions())
    afwMath.randomGaussianImage(randomImg, afwRandom)
    randomImg*= 20

    stats = afwMath.makeStatistics(randomImg, afwMath.VARIANCE)
    variance = stats.getValue(afwMath.VARIANCE)
    print variance

    img = modelImage.getImage()
    img += randomImg
    modelImage.getVariance().set(variance)

    exp = afwImage.ExposureF(modelImage, wcs)
    exp.setPsf(psf)

    ds9.mtv(exp, frame = 0)

    expList = measMult.ExposureListF()
    expList.append(exp)
    
    jiggeredLogShear = afwGeom.ellipses.LogShear(logShear[0]*1.1, logShear[1]*1.1, logShear[2]*1.1)
    #flux *= 1.1
    testModel = measMult.createSersicModel(flux, crVal, jiggeredLogShear)
    #testModel = model
    modelEvaluator = measMult.ModelEvaluator(testModel)
    modelEvaluator.setExposureList(expList)

    fitterPolicy = pexPolicy.Policy()
    fitterPolicy.add("checkGradient", True)

    fitter = measMult.MinuitAnalyticFitter(fitterPolicy)

    errors = [0.1, 0.1, 0.05, 0.1, 0.1, 0.1]
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
    fp = om.computeProjectionFootprint(psf, wcs)
    bbox = fp.getBBox()
    proj = om.makeProjection(psf, wcs, fp)    
    outputImage = afwImage.MaskedImageF(bbox.getWidth(), bbox.getHeight())
    outputImage.setXY0(bbox.getX0(), bbox.getY0())

    imageVector = proj.computeModelImage()
    varianceVector = numpy.zeros(fp.getNpix(), dtype=numpy.float32)
  
    measMult.expandImageF(fp, outputImage, imageVector, varianceVector)
    oe = afwImage.ExposureF(outputImage, wcs)
    ds9.mtv(oe, frame=1)



    
if __name__ == "__main__":
    applyFitter()
