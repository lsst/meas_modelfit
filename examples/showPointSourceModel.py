#!/usr/bin/env python

import lsst.meas.multifit as mf
import lsst.meas.algorithms
import lsst.afw.image
import lsst.afw.detection
import lsst.afw.geom
import numpy

from makeImageStack import makeImageStack

def main():

    factory = mf.PointSourceModelFactory()
    position = lsst.afw.geom.makePointD(0, 0)
    model = factory.makeModel(1.0,position)
    exposureList = makeImageStack(model, 1, position.getX(), position.getY())

    exposure = exposureList.front()
    exposure.writeFits("PointSourceProjection")
    
if __name__== "__main__":
    main()
