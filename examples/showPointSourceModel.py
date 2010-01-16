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
    position = lsst.afw.geom.makePointD(45, 45)
    model = factory.makeModel(1.0,position)
    exposureList = makeImageStack(model, 1, 45, 45)

    exposure = exposureList.front()
    exposure.writeFits("modelProjection")
    
if __name__== "__main__":
    main()
