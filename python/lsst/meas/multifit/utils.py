import lsst.afw.detection as afwDetection
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


def makeEllipseCore(src):    
    ixx= src.getIxx()
    ixy= src.getIxy()
    iyy= src.getIyy()

    return  afwGeom.ellipses.Quadrupole(ixx, iyy, ixy)

def checkEllipseCore(core, fp):
    axes = afwGeom.ellipses.Axes(core)
    if(numpy.isnan(axes.getA()) or numpy.isnan(axes.getB())):
        raise RuntimeError("Initial source measurements are NaN")

    maxR = 2*numpy.sqrt(fp.getArea())

    #print axes.getA(), axes.getB()
    if (axes.getA() <= 0.):
        axes.setA(1e-16)
    if (axes.getB() <= 0.):
        axes.setB(1e-16)
    if (axes.getA() > maxR or axes.getB() > maxR):
        raise RuntimeError("Initial source moments are unreasonably large")

    return axes

def makePoint(src):
    return afwGeom.Point2D(src.getXAstrom(), src.getYAstrom())


def checkSrcFlags(src):
    return True

def makeBitMask(mask, maskPlaneNames):
    bitmask=0
    for name in maskPlaneNames:
        bitmask |= mask.getPlaneBitMask(name)    
    return bitmask

def fitSource(exposure, src, bitmask, policy):
    if not checkSrcFlags(src):
        #print "Ignoring flagged source"
        return (None, None)

    #ftol = policy.get("ftol")
    #gtol = policy.get("gtol")
    #minStep = policy.get("minStep")
    #maxIter = policy.get("maxIter")
    #tau = policy.get("tau")
    #useSVD = policy.get("useSVD")
    
    fp = afwDetection.growFootprint(src.getFootprint(), policy.get("nGrowFp"))

    point = makePoint(src)
    try:
        core = makeEllipseCore(src)
        core = checkEllipseCore(core, fp)
        ellipse = afwGeom.ellipses.Ellipse(core, point)
        del core
    except Exception, e:
        print e
        ellipse = None

    optimizer = multifitLib.GaussNewtonOptimizer

    #print cutout
    #print bitmask
    #print fp
    #print ellipse
    #print point

    psDef = multifitLib.Definition.make(
           exposure, fp, point, 
           policy.getBool("isVariable"), 
           policy.getBool("isPositionActive"), 
           bitmask)
    psEval = multifitLib.Evaluator.make(psDef)

    try:
        psDistribution = optimizer.solve(psEval)
    except:
        psDistribution = None

    if not psDistribution:
        psInterpreter = None
    else:
        psInterpreter = multifitLib.UnifiedSimpleInterpreter.make(\
                psDistribution, psEval.getGrid())
    
    if not ellipse:
        return (psInterpreter, None)

    basis = loadBasis(policy.get("basisName"))
    sgDef = multifitLib.Definition.make(\
            exposure, fp, basis, ellipse,
            policy.get("isEllipticityActive"),
            policy.get("isRadiusActive"),
            policy.get("isPositionActive"),
            bitmask)
    sgEval = multifitLib.Evaluator.make(sgDef)
    try:
        sgDistribution = optimizer.solve(sgEval)
    except:
        sgDistribution = None

    if not sgDistribution:
        sgInterpreter = None
    else:
        sgInterpreter = multifitLib.UnifiedSimpleInterpreter.make(
            sgDistribution,
            sgEval.getGrid())

    return (psInterpreter, sgInterpreter) 


def processExposure(exposure, sources, policy):
    bitmask = makeBitmask(exposure.getMaskedImage().getMask(),
                          policy.getArra("maskPlaneName"))
    results = []   
    for s in sources:
        print "Fitting source", s.getId()
        results.append(fitSource(exposure, s, bitmask, policy))

    return results
