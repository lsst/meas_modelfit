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

import lsst.meas.multifit as measMult
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as geomEllipses
import lsst.afw.math as afwMath
import lsst.afw.detection as afwDet
import sys
import lsst.pex.policy as pexPol

import numpy
import lsst.afw.display.ds9 as ds9
import eups
import math

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
    print "variance", variance
    variance = 0.25
    modelImage.getVariance().set(variance)

    return modelImage

def main():
    numpy.set_printoptions(threshold=numpy.nan)
    i =0
    centroid = afwGeom.makePointD(45,45)
    psf = afwDet.createPsf("DoubleGaussian", 9, 9, 1.0)
    affine = afwGeom.AffineTransform()

    flux = 35.0
    axes = geomEllipses.Axes(30,25,0).transform(affine)
    print >> sys.stderr, axes

    logShear = geomEllipses.LogShear(axes)

    sersicIndex = 1.25
    pol = pexPol.Policy()
    cache = measMult.makeRobustSersicCache(pol)
    measMult.SersicMorphology.setSersicCache(cache)
    model = measMult.createSersicModel(flux, centroid, logShear, sersicIndex)
    #model = measMult.createPointSourceModel(flux, centroid)

    fp = model.computeProjectionFootprint(psf, affine)
    proj = model.makeProjection(psf, affine, fp)
    box = fp.getBBox()
    modelImage = afwImage.MaskedImageF(box.getWidth(), box.getHeight())
    modelImage.setXY0(box.getX0(), box.getY0())
    
    imageVector = proj.computeModelImage()
    varianceVector = numpy.zeros_like(imageVector)
    measMult.expandImageF(fp, modelImage, imageVector, varianceVector)

    #exp = makeModelExposure(model, psf, wcs)
    #expList = measMult.ExposureListF()
    #expList.append(exp)


    #eval = measMult.ModelEvaluator(model)
    #eval.setExposureList(expList)
    
    analyticDerivative = proj.computeNonlinearParameterDerivative()
    nonlinearParameters = model.getNonlinearParameters()
    
    #analyticDerivative = eval.computeNonlinearParameterDerivative()
    #nonlinearParameters = eval.getNonlinearParameters()

    eps = 2.2e-16

    
    columnVectors = []
    for n in range(model.getNonlinearParameterSize()):
        h = math.sqrt(eps) * nonlinearParameters[n]
        nonlinearParameters[n] += h        
        print  >> sys.stderr, "n: %g"%n
        print  >> sys.stderr, "h: %g"%h
        print >> sys.stderr, "params plus h: %s"%nonlinearParameters
        model.setNonlinearParameters(nonlinearParameters)
        print >> sys.stderr, model.getNonlinearParameters()

        plus = numpy.copy(proj.computeModelImage())

        nonlinearParameters[n] -= 2*h

        print >> sys.stderr, "params minus h: %s"%nonlinearParameters
        model.setNonlinearParameters(nonlinearParameters)
        minus = numpy.copy(proj.computeModelImage())

        partial = (plus - minus)
        partial = partial / (2*h)
        
        var = numpy.zeros_like(partial.base[0])
        fpBox = fp.getBBox()
        derivativeImage = afwImage.MaskedImageD(fpBox.getWidth(), fpBox.getHeight());
        derivativeImage.setXY0(fpBox.getLLC())
        measMult.expandImageD(fp, derivativeImage, numpy.array(analyticDerivative)[n, : ], var)
        ds9.mtv(derivativeImage, frame=n)

        columnVectors.append(partial)
        nonlinearParameters[n] += h
   
   
    numericDerivative = numpy.concatenate(columnVectors)
    diff = analyticDerivative / numericDerivative

    #print >> sys.stderr, "analytic:\n%s"% analyticDerivative
    #print >> sys.stderr, "numeric:\n%s"% numericDerivative
    #print >> sys.stderr, "diff:\n%s"%diff

    
if __name__== "__main__":
    main()
