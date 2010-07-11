#!/usr/bin/env python

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
