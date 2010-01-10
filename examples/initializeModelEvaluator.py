import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg
import lsst.meas.multifit as measMult
import numpy
import numpy.random

from makeImageStack import makeImageStack

def initializeModelEvaluator():

    psFactory = measMult.PointSourceModelFactory()
    psModel = psFactory.makeModel(1.0, afwGeom.makePointD(0,0))

    exposureList = makeImageStack(psModel, 15)
    modelEvaluator = measMult.ModelEvaluator(psModel, exposureList)

    print modelEvaluator.getNProjections()

if __name__ == "__main__":
    initializeModelEvaluator()
