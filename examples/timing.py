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

import lsst.meas.multifit as measMult
import lsst.afw.geom as afwGeom
import lsst.pex.policy as pexPolicy
import time

from makeImageStack import makeImageStack

def psTiming(maxIteration, depth):
    flux = 1.0
    centroid = afwGeom.makePointD(0,0)
    psModel = measMult.createPointSourceModel(flux, centroid)
    exposureList = makeImageStack(psModel, depth, centroid[0], centroid[1])
    t0 = time.time()
    modelEvaluator = measMult.ModelEvaluator(psModel)
    modelEvaluator.setExposureList(exposureList)
    t1 = time.time()
    print "Construction of ModelEvaluator: %0.3f ms"%((t1-t0)*1000.0)
    print "\tUsing %d pixels in %d exposures"% \
        (modelEvaluator.getNPixels(), modelEvaluator.getNProjections())

    t0 = time.time()
    fitterPolicy = pexPolicy.Policy()
    fitterPolicy.add("terminationType", "iteration")
    fitterPolicy.add("iterationMax", maxIteration)
    fitter=measMult.SingleLinearParameterFitter(fitterPolicy)
    t1 = time.time()
    print "Construction of Fitter: %0.3f ms"%((t1-t0)*1000.0)

    t0 = time.time()
    result = fitter.apply(modelEvaluator)
    t1 = time.time()
    nIterations = result.sdqaMetrics.get("nIterations")
    print "%d iterations of Fitter: %0.3f ms"%(nIterations, (t1-t0)*1000.0)

if __name__ == "__main__":
    psTiming(5, 30)
    print ""
    psTiming(10, 30)
    print ""
    psTiming(15, 30)
    print ""
    psTiming(20, 30)
    print ""
    psTiming(5, 50)
    print ""
    psTiming(20, 50)
    print ""
    psTiming(5, 200)
    print ""
    psTiming(20, 200)




