#!/usr/bin/env python

#
# LSST Data Management System
# Copyright 2008-2014 LSST Corporation.
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

from __future__ import print_function
from builtins import zip
import sys
import numpy

import lsst.meas.base
import lsst.meas.modelfit
import lsst.pex.logging

lsst.pex.logging.Debug("meas.modelfit.optimizer", 0)


def printResidualStatistics(image, msf, index):
    residuals = image.Factory(image, True)
    residuals *= -1
    msf.evaluate().addToImage(residuals)
    r = residuals.getArray()
    print("  %s: abs_diff_max=%f rms_diff=%f" % (index, numpy.abs(r).max(), ((r**2).mean())**0.5))


def fitPsfImage(image, fitters):
    maskedImage = lsst.afw.image.MaskedImageD(image)
    result = lsst.meas.base.SdssShapeAlgorithm.Result()
    try:
        lsst.meas.base.SdssShapeAlgorithm.apply(
            maskedImage,
            lsst.afw.detection.Footprint(image.getBBox(lsst.afw.image.PARENT)),
            lsst.afw.geom.Point2D(0.0, 0.0),
            result
        )
    except lsst.meas.base.MeasurementError:
        pass
    if (result.x**2 + result.y**2) > 0.25:
        raise ValueError("Bad centroid: %s" % result.getCentroid())
    msf = fitters[0].apply(image, result.getShape())
    n = 0
    for previousFitter, nextFitter in zip(fitters[:-1], fitters[1:]):
        printResidualStatistics(image, msf, n)
        msf = nextFitter.adapt(msf, previousFitter.getModel())
        msf = nextFitter.apply(image, msf)
        n += 1
    printResidualStatistics(image, msf, n)

ctrls = [
    lsst.meas.modelfit.GeneralPsfFitterControl(),
    lsst.meas.modelfit.GeneralPsfFitterControl(),
]
ctrls[0].primary.ellipticityPriorSigma = 0.3
ctrls[0].primary.radiusPriorSigma = 0.5
ctrls[0].wings.ellipticityPriorSigma = 0.3
ctrls[0].wings.radiusPriorSigma = 0.5
ctrls[1].inner.order = 0
ctrls[1].primary.order = 4
ctrls[1].wings.order = 4
ctrls[1].outer.order = 0
ctrls[1].inner.ellipticityPriorSigma = 0.3
ctrls[1].inner.radiusPriorSigma = 0.5
ctrls[1].inner.positionPriorSigma = 0.1
ctrls[1].primary.ellipticityPriorSigma = 0.3
ctrls[1].primary.radiusPriorSigma = 0.5
ctrls[1].primary.positionPriorSigma = 0.1
ctrls[1].wings.ellipticityPriorSigma = 0.3
ctrls[1].wings.radiusPriorSigma = 0.5
ctrls[1].wings.positionPriorSigma = 0.1
ctrls[1].outer.ellipticityPriorSigma = 0.3
ctrls[1].outer.radiusPriorSigma = 0.5
ctrls[1].outer.positionPriorSigma = 0.1


def main(filenames):
    fitters = [lsst.meas.modelfit.GeneralPsfFitter(ctrl) for ctrl in ctrls]
    for filename in filenames:
        image = lsst.afw.image.ImageD(filename)
        image.getArray()[:, :] /= image.getArray()[:, :].max()
        print(filename)
        fitPsfImage(image, fitters)

if __name__ == "__main__":
    main(sys.argv[1:])
