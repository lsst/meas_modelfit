import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses
from . import multifitLib
import numpy
import os
import eups

def loadBasis(name):
    productDir = eups.productDir("meas_multifit")
    path = os.path.join(productDir, "data", "%s.boost" % name)
    return multifitLib.CompoundShapeletModelBasis.load(path)


def makeEllipse(src):    
    ixx= src.getIxx()
    ixy= src.getIxy()
    iyy= src.getIyy()
    ixx = ((ixx < 0 or numpy.isnan(ixx)) and [0.0] or [ixx])[0]
    ixy = ((ixy < 0 or numpy.isnan(ixy)) and [0.0] or [ixy])[0]
    ixx = ((iyy < 0 or numpy.isnan(iyy)) and [0.0] or [iyy])[0]

    quad = afwGeom.ellipses.Quadrupole(
        src.getIxx(), 
        src.getIyy(), 
        src.getIxy()
    )
    point = makePoint(src)
    return afwGeom.ellipses.Ellipse(quad, point)

def makePoint(src):
    return afwGeom.Point2D(src.getXAstrom(), src.getYAstrom())


def checkSrcFlags(src):
    return True

def makeBitMask(mask, maskPlaneNames):
    bitmask=0
    for name in maskPlaneNames:
        bitmask |= mask.getPlaneBitMask(name)    


def fitSource(cutout, src, policy):
    optimizer = multifitLib.GaussNewtonOptimizer()

    bitmask = makeBitMask(cutout.getMaskedImage().getMask(), 
            policy.getArray("maskPlaneName"))

    ellipse = makeEllipse(src)
    fp = src.getFootprint()
    psDef = multifitLib.Definition.make(
           cutout, fp, ellipse.getCenter(), 
           policy.get("isVariable"), 
           policy.get("isPositionActive"), 
           bitmask)
    psEval = multifitLib.Evaluator(psDef)
    psDistribution = optimizer.solve(psEval)
    psInterpreter = multifitLib.SimpleUnifiedInterpreter.make(
            psDistribution,
            psEval.getGrid())
    
    basis = multifitLib.loadBasis(policy.get("basisName"))
    sgDef = multifitLib.Definition.make(
            cutout, fp, basis, ellipse,
            policy.get("isEllipticityActive"),
            policy.get("isRadiusActive"),
            policy.get("isPositionActive"),
            bitmask)
    sgEval = multifitLib.Evaluator(sgDef)
    sgDistribution = optimizer.solve(sgEval)
    sgInterpreter = multifitLib.SimpleUnifiedInterpreter.make(
            sgDistribution,
            sgEval.getGrid())

    return (psInterpreter, sgInterpreter) 

