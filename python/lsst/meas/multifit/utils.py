import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses
from . import multifitLib
import numpy
import os
import eups
import math

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
    axes = afwGeom.ellipses.Axes(quad)
    maxR = 2*numpy.sqrt(src.getFootprint().getArea())

    print axes.getA(), axes.getB()
    if (axes.getA() == 0.):
        axes.setA(1e-16)
    if (axes.getB() == 0.):
        axes.setB(1e-16)
    if (axes.getA() > maxR or axes.getB() > maxR or \
            numpy.isnan(axes.getA()) or numpy.isnan(axes.getB())):
        raise RuntimeError("Initial Source moments are unreasonably large")

    point = makePoint(src)
    return afwGeom.ellipses.Ellipse(axes, point)

def makePoint(src):
    return afwGeom.Point2D(src.getXAstrom(), src.getYAstrom())


def checkSrcFlags(src):
    return True

def makeBitMask(mask, maskPlaneNames):
    bitmask=0
    for name in maskPlaneNames:
        bitmask |= mask.getPlaneBitMask(name)    
    return bitmask

def fitSource(cutout, src, policy):
    optimizer = multifitLib.GaussNewtonOptimizer()

    bitmask = makeBitMask(cutout.getMaskedImage().getMask(), 
            policy.getArray("maskPlaneName"))

    #ftol = policy.get("ftol")
    #gtol = policy.get("gtol")
    #minStep = policy.get("minStep")
    #maxIter = policy.get("maxIter")
    #maxIter = policy.get("tau")
    #maxIter = policy.get("useSVD")
    
    try:
        ellipse = makeEllipse(src)
        point = ellipse.getCenter()
    except Exception, e:
        print e
        ellipse = None
        point = makePoint(src)        


    fp = src.getFootprint()

    print cutout
    print bitmask
    print fp
    print ellipse
    print point

    psDef = multifitLib.Definition.make(
           cutout, fp, point, 
           policy.getBool("isVariable"), 
           policy.getBool("isPositionActive"), 
           bitmask)
    psEval = multifitLib.Evaluator.make(psDef)
    psDistribution = optimizer.solve(psEval)

    if not psDistribution:
        psInterpreter = None
    else:
        psInterpreter = multifitLib.UnifiedSimpleInterpreter.make(\
                psDistribution, psEval.getGrid())
    
    if not ellipse:
        return (psInterpreter, None)

    basis = loadBasis(policy.get("basisName"))
    sgDef = multifitLib.Definition.make(\
            cutout, fp, basis, ellipse,
            policy.get("isEllipticityActive"),
            policy.get("isRadiusActive"),
            policy.get("isPositionActive"),
            bitmask)
    sgEval = multifitLib.Evaluator.make(sgDef)
    sgDistribution = optimizer.solve(sgEval)

    if not sgDistribution:
        sgInterpreter = None
    else:
        print sgEval.getParameterSize()
        print sgEval.getCoefficientSize()
        sgInterpreter = multifitLib.UnifiedSimpleInterpreter.make(
            sgDistribution,
            sgEval.getGrid())

    return (psInterpreter, sgInterpreter) 

