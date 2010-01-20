import sys
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as geomEllipse
import lsst.afw.image as afwImg
import lsst.afw.detection as afwDet
import lsst.meas.multifit as measMult

def footprintFromEllipse(x, y, a, b, theta, filename="footprint.fits"):
    axes = geomEllipse.Axes(a, b, theta)
    ellipse = geomEllipse.AxesEllipse(axes, afwGeom.makePointD(x, y))
    fp = measMult.makeFootprint(ellipse)
    bbox = fp.getBBox()
    mask = afwImg.MaskU(bbox.getWidth(), bbox.getHeight())
    mask.setXY0(bbox.getX0(), bbox.getY0())
    
    afwDet.setMaskFromFootprint(mask, fp, 1)
    mask.writeFits(filename)
    return fp

if __name__=='__main__':
    argc = len(sys.argv)

    if(argc == 6):
        footprintFromEllipse(
            sys.argv[1], sys.argv[2], 
            sys.argv[3], sys.argv[4], sys.argv[5]
        )
    elif (len(sys.argv) == 7):
        footprintFromEllipse(
            sys.argv[1], sys.argv[2], 
            sys.argv[3], sys.argv[4], sys.argv[5],
            sys.argv[6]
        )
    else:
        print "wrong number of args. x, y, a, b, theta [, filename]"

