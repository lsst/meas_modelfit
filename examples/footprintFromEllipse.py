import sys
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as geomEllipse
import lsst.afw.image as afwImg
import lsst.afw.detection as afwDet
import lsst.meas.multifit as measMult
import lsst.afw.display.ds9 as ds9

def footprintFromEllipse(x, y, a, b, theta):
    axes = geomEllipse.Axes(float(a), float(b), float(theta))
    point = afwGeom.makePointD(float(x), float(y))
    ellipse = geomEllipse.AxesEllipse(axes, point)
    fp = measMult.makeFootprint(ellipse)

    return fp
