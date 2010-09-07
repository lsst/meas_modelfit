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


def applyFitter():
    #exp = afwImage.ExposureF("c00.fits")
    #wcs = exp.getWcs()
    crVal = afwGeom.makePointD(0,0)
    crPix = afwGeom.makePointD(1,1)
    wcs = afwImage.createWcs(crVal, crPix, 1., 0., 0., 1.)

    axes = afwGeom.ellipses.Axes(25,30,0)
    affine = wcs.linearizeAt(crVal)

    transformedAxes = axes.transform(affine)
    logShear = afwGeom.ellipses.LogShear(transformedAxes)
    sersicIndex = 1.5
    flux = 1.0

    model = measMult.createSersicModel(flux, crVal, logShear, sersicIndex)

    psf = afwDet.createPsf("DoubleGaussian", 7, 7, 1.0)

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
    afwMath.randomUniformImage(randomImg, afwRandom)
    
    img = modelImage.getImage()
    img += randomImg
    modelImage.getVariance().set(0.25)

    exp = afwImage.ExposureF(modelImage, wcs)
    exp.setPsf(psf)

    ds9.mtv(exp, frame = 0)

    expList = measMult.ExposureListF()
    expList.append(exp)

    modelEvaluator = measMult.ModelEvaluator(model)
    modelEvaluator.setExposureList(expList)
    
    fitterPolicy = pexPolicy.Policy()
    fitterPolicy.add("terminationType", "iteration")
    fitterPolicy.add("terminationType", "dChisq")
    fitterPolicy.add("iterationMax", 5)
    fitterPolicy.add("dChisqThreshold", 0.0001)

    fitter = measMult.SingleLinearParameterFitter(fitterPolicy)
    result = fitter.apply(modelEvaluator)

    print "nIterations: %d"%result.sdqaMetrics.get("nIterations")
    print "chisq: %d"%result.chisq
    print "dChisq: %d"%result.dChisq
    print modelEvaluator.getLinearParameters()
    print modelEvaluator.getNonlinearParameters()
   
    om = modelEvaluator.getModel().clone()
    #fp = om.computeProjectionFootprint(psf, wcs)
    #bbox = fp.getBBox()
    proj = om.makeProjection(psf, wcs, fp)    
    outputImage = afwImage.MaskedImageF(grownBbox.getWidth(), grownBbox.getHeight())
    outputImage.setXY0(grownBbox.getX0(), grownBbox.getY0())

    imageVector = proj.computeModelImage()
    varianceVector = numpy.zeros(nPix, dtype=numpy.float32)
  
    subImage = outputImage.Factory(outputImage, bbox)
    measMult.expandImageF(fp, outputImage, imageVector, varianceVector)
    oe = afwImage.ExposureF(outputImage, wcs)
    ds9.mtv(oe, frame=1)

    
if __name__ == "__main__":
    applyFitter()
