#!/usr/bin/env python

import lsst.meas.multifit as mf
import lsst.meas.algorithms
import lsst.afw.image
import lsst.afw.detection
import lsst.afw.geom
import numpy

def main():
    crVal = lsst.afw.image.PointD(45,45)
    crPix = lsst.afw.image.PointD(0,0)    
    wcs = lsst.afw.image.createWcs(crVal, crPix, 0.0001, 0, 0, 0.0001)
    psf = lsst.meas.algorithms.createPSF("DoubleGaussian", 19,19, 2)
    factory = mf.PointSourceModelFactory()


    
    position = lsst.afw.geom.makePointD(45, 45)
    model = factory.makeModel(1.0,position)
    footprint = model.computeProjectionFootprint(psf, wcs)
    bbox = footprint.getBBox()
    proj = model.makeProjection(psf,wcs,footprint)
    modelImage = proj.computeModelImage()
    varianceVector = numpy.ones(footprint.getNpix(), dtype=numpy.float32)
    varianceVector *= 0.25
    exp = mf.CharacterizedExposureD(bbox.getWidth(), bbox.getHeight(), wcs, psf)
    maskedImage = exp.getMaskedImage()
    maskedImage.setXY0(bbox.getX0(), bbox.getY0())
    mf.expandImageD(footprint, maskedImage, modelImage, varianceVector)
    lsst.afw.detection.setMaskFromFootprint(maskedImage.getMask(), footprint, 1)
    maskedImage.writeFits("modelProjection")

    print proj.computeModelImage()
    print proj.computeLinearParameterDerivative()
    print proj.computeNonlinearParameterDerivative()

if __name__== "__main__":
    main()
