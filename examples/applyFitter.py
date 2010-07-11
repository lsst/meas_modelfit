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

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg
import lsst.meas.multifit as measMult
import lsst.pex.policy as pexPolicy
import numpy
import numpy.random

from makeImageStack import makeImageStack

def applyFitter():
    centroid = afwGeom.makePointD(45,45)
    flux = 1.0
    psModel = measMult.createPointSourceModel(flux, centroid)

    exposureList = makeImageStack(psModel, 15, 45, 45)
    modelEvaluator = measMult.ModelEvaluator(psModel, exposureList)
    
    fitterPolicy = pexPolicy.Policy()
    fitterPolicy.add("terminationType", "iteration")
    fitterPolicy.add("terminationType", "dChisq")
    fitterPolicy.add("iterationMax", 5)
    fitterPolicy.add("dChisqThreshold", 0.0001)

    fitter = measMult.SingleLinearParameterFitter(fitterPolicy)
    result = fitter.apply(modelEvaluator)

    print "nIterations: %d"%result.sdqaMetrics.get("nIterations")
    print "chisq: %d"%result.chisq
    print "dChisq: %d"%result.dChisq

if __name__ == "__main__":
    applyFitter()
