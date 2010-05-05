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
