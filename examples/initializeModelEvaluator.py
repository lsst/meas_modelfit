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
    print "Number of exposures in ModelEvaluator: %d"%modelEvaluator.getNProjections()
    print "Number of total pixels in ModelEvaluator: %d"%modelEvaluator.getNPixels()

    print "ModelEvaluator image vector:"
    print modelEvaluator.getImageVector()

    print "Model Evaluator ModelImage:"    
    print modelEvaluator.computeModelImage()

    print "Model Evaluator LinearParameterDerivative:"    
    print modelEvaluator.computeLinearParameterDerivative()

    print "Model Evaluator NoninearParameterDerivative:"    
    print modelEvaluator.computeNonlinearParameterDerivative()


if __name__ == "__main__":
    initializeModelEvaluator()
