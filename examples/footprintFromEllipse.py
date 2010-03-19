import sys
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as geomEllipse
import lsst.afw.image as afwImg
import lsst.afw.detection as afwDet
import lsst.meas.multifit as measMult

def footprintFromEllipse(x, y, a, b, theta, filename="footprint.fits"):
    axes = geomEllipse.Axes(float(a), float(b), float(theta))
    point = afwGeom.makePointD(float(x), float(y))
    ellipse = geomEllipse.AxesEllipse(axes, point)
    fp = measMult.makeFootprint(ellipse)
    bbox = fp.getBBox()
    mask = afwImg.MaskU(bbox.getWidth(), bbox.getHeight())
    mask.setXY0(bbox.getLLC())
    
    afwDet.setMaskFromFootprint(mask, fp, 1)
    mask.writeFits(filename)

    return fp

if __name__=='__main__':
    argc = len(sys.argv)

    if(argc == 6 or argc == 7):
        footprintFromEllipse(*sys.argv[1:])
    else:
        print "wrong number of args. x, y, a, b, theta [, filename]"

