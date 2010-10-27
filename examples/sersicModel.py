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
import lsst.afw.geom.ellipses
import lsst.afw.math as afwMath
import lsst.afw.detection as afwDet

import numpy
import lsst.afw.display.ds9 as ds9
import eups

def makeModelExposure(model, psf, wcs, noiseFactor=0):
    fp = model.computeProjectionFootprint(psf, wcs)

    box = fp.getBBox()
    grow = afwImage.BBox(box.getLLC(), int(box.getWidth()*5), int(box.getHeight()*5))
    grow.shift(-int(box.getWidth()*2), -int(box.getHeight()*2))
    growFp = afwDet.Footprint(grow)
    nPix = growFp.getNpix()

    proj = model.makeProjection(psf, wcs, growFp)
    modelImage = afwImage.MaskedImageF(grow.getWidth(), grow.getHeight())
    modelImage.setXY0(grow.getX0(), grow.getY0())
    
    imageVector = proj.computeModelImage()
    varianceVector = numpy.zeros(nPix, dtype=numpy.float32)

    measMult.expandImageF(growFp, modelImage, imageVector, varianceVector)

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

def main():
    i =0
    crVal = afwGeom.makePointD(45,45)
    crPix = afwGeom.makePointD(1,1)
    wcs = afwImage.createWcs(crVal, crPix, 0.0001, 0., 0., 0.0001)

    psf = afwDet.createPsf("DoubleGaussian", 9, 9, 1.0)
    affine = wcs.linearizePixelToSky(crVal)

    flux = 5.0
    centroid = crVal
    axes = lsst.afw.geom.ellipses.Axes(30,25,0).transform(affine)
    print axes

    logShear = lsst.afw.geom.ellipses.LogShear(axes)

    sersicIndex = 1.25

    try:
        root = os.path.join(eups.productDir("multifitData"), "cache")
        path = os.path.join(root, "sersicCache.boost")
        cache = measMult.SersicCache.load(path)
    except:
        pol = pexPol.Policy()
        cache = measMult.SersicCache.make(pol)        
    measMult.SersicMorphology.setSersicCache(cache)

    model = measMult.createSersicModel(flux, centroid, logShear, sersicIndex)
    exp = makeModelExposure(model, psf, wcs)

    ds9.mtv(exp, frame=i)
    i+=1
    del model
    del exp

    sersicIndex = 4.0
    model = measMult.createSersicModel(flux, centroid, logShear, sersicIndex)
    exp = makeModelExposure(model, psf, wcs)
    ds9.mtv(exp, frame=i)
    i+=1
    del model
    del exp

    axes = lsst.afw.geom.ellipses.Axes(50,10,30).transform(affine)
    print axes
    logShear = lsst.afw.geom.ellipses.LogShear(axes)
    model = measMult.createSersicModel(flux, centroid, logShear, sersicIndex)
    exp = makeModelExposure(model, psf, wcs)
    ds9.mtv(exp, frame=i)

if __name__== "__main__":
    main()
