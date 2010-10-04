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


def makeModelExposure(model, psf, affine, noiseFactor=0):
    fp = model.computeProjectionFootprint(psf, affine)
    nPix = fp.getNpix()
    box = fp.getBBox()
    proj = model.makeProjection(psf, affine, fp)
    modelImage = afwImage.MaskedImageF(box.getWidth(), box.getHeight())
    modelImage.setXY0(box.getX0(), box.getY0())
    
    imageVector = proj.computeModelImage()
    varianceVector = numpy.zeros_like(imageVector)

    measMult.expandImageF(fp, modelImage, imageVector, varianceVector)

    afwRandom = afwMath.Random()
    randomImg = afwImage.ImageF(modelImage.getDimensions())
    afwMath.randomGaussianImage(randomImg, afwRandom)
    randomImg *= noiseFactor
    img = modelImage.getImage()
    img += randomImg

    stats = afwMath.makeStatistics(modelImage, afwMath.VARIANCE)
    variance = stats.getValue(afwMath.VARIANCE)
    modelImage.getVariance().set(variance)
    
    return modelImage

def applyFitter():
    #exp = afwImage.ExposureF("c00.fits")
    #wcs = exp.getWcs()
    frameId = 0
    psf = afwDet.createPsf("DoubleGaussian", 9, 9, 1.0)
    affine = afwGeom.AffineTransform()
   
    pixel = afwGeom.makePointD(45, 45)
    axes = afwGeom.ellipses.Axes(25,30,0)
    sersicIndex = 1.5
    flux = 1.0

    model = measMult.createPointSourceModel(flux, pixel)
    errors = [0.1, 1e-5, 1e-5]

    exp = makeModelExposure(model, psf, affine, 1.0)
    ds9.mtv(exp, frame=frameId)
    frameId +=1

    modelEvaluator = measMult.ModelEvaluator(model, affine)
    modelEvaluator.setData(exp, psf, affine)

    fitterPolicy = pexPolicy.Policy()
    fitterPolicy.set("checkGradient", True)
    fitter = measMult.MinuitAnalyticFitter(fitterPolicy)

    result = fitter.apply(modelEvaluator, errors)

    print "nPix", modelEvaluator.getNPixels()
    print "nIterations", result.nIterations
    print "chisq", result.chisq

    om = modelEvaluator.getModel().clone()
    outExp = makeModelExposure(om, psf, affine)

    ds9.mtv(outExp, frame=frameId)



    
if __name__ == "__main__":
    applyFitter()
