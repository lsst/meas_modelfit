import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg
import lsst.meas.multifit as measMult
import numpy
import numpy.random

from makeImageStack import makeImageStack

def initializeModelEvaluator():
    psFactory = measMult.PointSourceModelFactory()
    psModel = psFactory.makeModel(1.0, afwGeom.makePointD(45,45))

    exposureList = makeImageStack(psModel, 15, 45, 45)
    modelEvaluator = measMult.ModelEvaluator(psModel, exposureList)

    numpy.set_printoptions(threshold=100000)
    print "ModelEvaluator nProjections: %d"%modelEvaluator.getNProjections()
    print "ModelEvaluator nPixels: %d"%modelEvaluator.getNPixels()

    print "ModelEvaluator image vector: %s"%modelEvaluator.getImageVector()
    print "ModelEvaluator variance vector: %s"%modelEvaluator.getVarianceVector()  
    print "ModelEvaluator ModelImage: %s"%modelEvaluator.computeModelImage()
    print "ModelEvaluator LinearParameterDerivative: %s"%modelEvaluator.computeLinearParameterDerivative()
    print "ModelEvaluator NoninearParameterDerivative: %s"%modelEvaluator.computeNonlinearParameterDerivative()


if __name__ == "__main__":
    initializeModelEvaluator()
