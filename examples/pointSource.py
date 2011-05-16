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

import lsst.afw.image
import lsst.afw.detection
import lsst.afw.image
import lsst.afw.geom.ellipses
import lsst.meas.multifit
import numpy
from matplotlib import pyplot

def makeBasis():
    psf = lsst.afw.detection.createPsf("DoubleGaussian", 19, 19, 2.0, 1.0)
    localPsf = psf.getLocalPsf(lsst.afw.geom.Point2D(0.0, 0.0))
    basis = lsst.meas.multifit.loadBasis("ed+00:0000")
    convolvedBasis = basis.convolve(localPsf)
    return convolvedBasis

basis = makeBasis()
bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-10, -10), lsst.afw.geom.Point2I(10, 10))
footprint = lsst.afw.detection.Footprint(bbox)
envelope = lsst.afw.geom.Box2D(bbox)
extent = (envelope.getMinX(), envelope.getMaxX(), envelope.getMinY(), envelope.getMaxY())
core = lsst.meas.multifit.EllipseCore(5.0, 0.0, 1E-16)
ellipse = lsst.afw.geom.ellipses.Ellipse(core, lsst.afw.geom.Point2D(0.0, 0.0))
array = numpy.zeros((footprint.getArea(), basis.getSize()), dtype=float)
images = array.reshape(bbox.getHeight(), bbox.getWidth(), basis.getSize())
basis.evaluate(array, footprint, ellipse)
