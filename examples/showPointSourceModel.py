#!/usr/bin/env python

import lsst.meas.multifit as mf
import lsst.meas.algorithms
import lsst.afw.image
import lsst.afw.geom

def main():
    wcs = lsst.afw.image.Wcs()
    psf = lsst.meas.algorithms.createPSF("DoubleGaussian", 9,9, 3)
    factory = mf.PointSourceModelFactory()
    position = lsst.afw.geom.PointD()
    model = factory.makeModel(1.0,position)
    proj = model.makeProjection(psf,wcs,model.computeProjectionFootprint(psf,wcs))
    proj.computeModelImage()
    proj.computeLinearParameterDerivative()
    proj.computeNonlinearParameterDerivative()
    #proj.computeWcsParameterDerivative()
    proj.computePsfParameterDerivative()

if __name__== "__main__":
    main()
