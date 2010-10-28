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
import lsst.pex.policy as pexPol
import sys, os
import eups

def makeModelExposure(model, psf, affine, noiseFactor=0):
    fp = model.computeProjectionFootprint(psf, affine)

    box = fp.getBBox()

    del fp
    bigBox = afwImage.BBox(box.getLLC(), box.getWidth()+120, box.getHeight()+120);
    bigBox.shift(-60, -60);
    bigFP = afwDet.Footprint(bigBox);

    del box
    nPix = bigFP.getNpix()

    
    print >> sys.stderr, "\tmaking pojection"
    proj = model.makeProjection(psf, affine, bigFP)
    modelImage = afwImage.MaskedImageF(bigBox.getWidth(), bigBox.getHeight())
    modelImage.setXY0(bigBox.getX0(), bigBox.getY0())
    
    print >> sys.stderr, "\tmaking buffers"
    imageVector = proj.computeModelImage()
    varianceVector = numpy.zeros_like(imageVector)

    print >> sys.stderr, "\texpanding data"
    measMult.expandImageF(bigFP, modelImage, imageVector, varianceVector)


    print >> sys.stderr, "\tadding noise"
    afwRandom = afwMath.Random()
    print >> sys.stderr, "\t\tallocating random image"
    randomImg = afwImage.ImageF(modelImage.getDimensions())

    print >> sys.stderr, "\t\tcompute random image"
    afwMath.randomGaussianImage(randomImg, afwRandom)


    print >> sys.stderr, "\t\tscale by noise"
    randomImg *= noiseFactor


    print >> sys.stderr, "\t\tadd random to modelImage"
    img = modelImage.getImage()
    img += randomImg
    del randomImg

    print >> sys.stderr, "\tcomputing variance"
    stats = afwMath.makeStatistics(modelImage, afwMath.VARIANCE)
    variance = stats.getValue(afwMath.VARIANCE)

    del stats
    modelImage.getVariance().set(variance)
    
    return modelImage

def applyFitter():
    frameId = 0
    psf = afwDet.createPsf("DoubleGaussian", 9, 9, 1.0)
    affine = afwGeom.AffineTransform()
   
    pixel = afwGeom.makePointD(45, 45)
    axes = afwGeom.ellipses.Axes(25,30,0)
    sersicIndex = 1.5
    flux = 1.0
    

    try:
        print>>sys.stderr, "loading cache"
        root = os.path.join(eups.productDir("meas_multifitData"), "cache")
        path = os.path.join(root, "sersicCache.boost")
        cache = measMult.SersicCache.load(path)
    except:
        print >> sys.stderr, "...failed. making cache"
        pol = pexPol.Policy()
        cache = measMult.SersicCache.make(pol)        
    measMult.SersicMorphology.setSersicCache(cache)

    print >> sys.stderr, "making model"
    #model = measMult.createPointSourceModel(flux, pixel)
    #errors = [0.1, 1e-5, 1e-5]
    model = measMult.createSersicModel(flux, pixel, axes, 1.0)
    errors = [0.1, 1e-5, 1e-5, 1e-8, 1e-8, 1e-8, 1e-5]


    print >> sys.stderr, "making image"
    exp = makeModelExposure(model, psf, affine, 1.0)
    ds9.mtv(exp, frame=frameId)
    frameId +=1


    print >> sys.stderr, "making evaluator"
    modelEvaluator = measMult.ModelEvaluator(model)
    modelEvaluator.setData(exp, psf, affine)

    fitterPolicy = pexPolicy.Policy()
    fitter = measMult.MinuitNumericFitter(fitterPolicy)

    result = fitter.apply(modelEvaluator, errors)

    print "nPix", modelEvaluator.getNPixels()
    print "nIterations", result.nIterations
    print "chisq", result.chisq
    print "chisq/dog", result.chisq/(modelEvaluator.getNPixels()-7)

    om = modelEvaluator.getModel().clone()
    outExp = makeModelExposure(om, psf, affine)

    ds9.mtv(outExp, frame=frameId)



    
if __name__ == "__main__":
    applyFitter()
