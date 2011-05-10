import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses
import lsst.meas.multifit as measMult
import lsst.meas.multifit.definition
import lsst.meas.multifit.grid
import numpy

def makeEllipse(src):    
    ixx= src.getIxx()
    ixy= src.getIxy()
    iyy= src.getIyy()
    ixx = ((ixx < 0 or numpy.isnan(ixx)) and [0.0] or [ixx])[0]
    ixy = ((ixy < 0 or numpy.isnan(ixy)) and [0.0] or [ixy])[0]
    ixx = ((iyy < 0 or numpy.isnan(iy)) and [0.0] or [iyy])[0]

    quad = afwGeom.ellipses.Quadrupole(
        src.getIxx(), 
        src.getIyy(), 
        src.getIxy()
    )
    point = makePoint(src)
    return afwGeom.ellipses(quad, point)

def makePoint(src):
    return afwGeom.Point2D(src.getXAstrom(), src.getYAstrom())


def checkSrcFlags(src):
    return True

def makeBitMask(mask, maskPlaneNames):
    bitmask=0
    for name in maskPlaneNames:
        bitmask |= mask.getPlaneBitMask(name)    


def fitSource(cutout, fp, src, policy):
    optimizer = measMult.GaussNewtonOptimizer()

    bitmask = makeBitMask(cutout.getMaskedImage().getMask(), 
            policy.getArray("maskPlaneName"))

    ellipse = makeEllipse(src)

    definition = measMult.Definition.make(
           cutout, fp, src.getCenter(), 
           policy.get("isVariable"), 
           policy.get("isPositionActive"), 
           bitmask)
    evaluator = measMult.Evaluator(definition)
    distribution = optimizer.solve(evaluator)

    basis = measMult.loadBasis(policy.get("basisName"))
    definition = measMult.Definition.make(
            cutout, fp, basis, ellipse,
            policy.get("isEllipticityActive"),
            policy.get("isRadiusActive"),
            policy.get("isPositionActive"),
            bitmask)
    evaluator = measMult.Evaluator(definition)
    distribution = optimizer.solve(evaluator)

    

