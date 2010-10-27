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

import sys
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.detection as afwDet
import lsst.meas.algorithms as measAlg
import lsst.meas.multifit as measMult
import numpy
import numpy.random
import lsst.afw.display.ds9 as ds9

def makeImageStack(model, depth):
    psf = afwDet.createPsf("DoubleGaussian", 7, 7, 1.0)
    identity = afwGeom.AffineTransform()
    fp = model.computeProjectionFootprint(psf, identity)
    projection = model.makeProjection(psf, identity, fp)

    bbox = fp.getBBox()
    nPix = fp.getNpix();

    modelImage = afwImage.MaskedImageF(bbox.getWidth(), bbox.getHeight())
    modelImage.setXY0(bbox.getX0(), bbox.getY0())

    imageVector = projection.computeModelImage()
    varianceVector = numpy.zeros(nPix, dtype=numpy.float32)

    measMult.expandImageF(fp, modelImage, imageVector, varianceVector)

    sigma = 0.5
    sigmaSq = sigma**2

    miList = []
    psfList = []
    transformList = []

    afwRandom = afwMath.Random()

    for i in xrange(depth):
        #create an exposure whose size matches footprint size
        mi = afwImage.MaskedImageF( 
            bbox.getWidth(), bbox.getHeight()
        )
        mi.setXY0(bbox.getX0(), bbox.getY0())

        img = mi.getImage()
        afwMath.randomGaussianImage(img, afwRandom)

        img += modelImage.getImage()
        del img

        var = mi.getVariance()
        var.set(sigmaSq)
        del var

        miList.append(mi)
        psfList.append(psf)
        transformList.append(identity)
    return  (miList, PsfList, transformList)

def displayImageStack(depth):
    flux = 1.0
    position = afwGeom.makePointD(45,45)
    psModel = measmultifit.createPointSourceModel(flux, position)
    (miList, psfList, transformList) = makeImageStack(psModel, depth)

    for i, mi in enumerate(miList):
        ds9.mtv(mi, title=str(i))

if __name__ == "__main__":
    argc = len(sys.argv)
    if argc == 1:
        displayImageStack(1)
    else if argc == 2:
        displayImageStack(int(sys.argv[1]))
    else:      
        usage =
        """
        Usage:
            makeImageStack [<depth> [<output prefix>]]

        <depth> - number of images to produce [1, ) defaults to 1
        """

