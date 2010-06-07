#!/usr/bin/env python

import lsst.meas.multifit as mf
import lsst.meas.algorithms
import lsst.afw.image
import lsst.afw.detection
import lsst.afw.geom
import numpy
import lsst.afw.display.ds9 as ds9
from makeImageStack import makeImageStack

def main():
    flux = 1.0
    position = lsst.afw.geom.makePointD(0, 0)
    model = mf.createPointSourceModel(flux,position)
    exposureList = makeImageStack(model, 1, position.getX(), position.getY())

    exposure = exposureList.front()    
    mi = exposure.getMaskedImage()
    ds9.mtv(mi, frame=0, wcs=exposure.getWcs())
    del mi
    del exposure
    del exposureList
    del model
    del position
    
if __name__== "__main__":
    main()
