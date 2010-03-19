import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as geomEllipses
import lsst.afw.image as afwImg
import lsst.afw.detection as afwDet
import lsst.meas.multifit as measMult
import numpy


def printFp(fp, prefix):
    bbox = fp.getBBox()
    img = afwImg.MaskedImageD(bbox.getWidth(), bbox.getHeight())
    img.setXY0(bbox.getLLC())
    data=numpy.arange(fp.getNpix())
    measMult.expandImageD(fp, img, data, data)
    filename = prefix + "_%d"%fp.getNpix()
    img.getImage().writeFits(filename+"_img.fits")
    afwDet.setMaskFromFootprint(img.getMask(), fp, 1)
    img.getMask().writeFits(filename+"_msk.fits")
    
    del img
    del bbox
    del data

def main():
    bbox = afwImg.BBox(afwImg.PointI(0,0), 8, 8)
    fp = afwDet.Footprint(bbox)
    printFp(fp, "boxFp")

    del bbox
    del fp

    core = geomEllipses.Axes(23, 15, 1.30)
    point = afwGeom.makePointD(3.5,3.5)
    ellipse = geomEllipses.AxesEllipse(core, point)
    fp = measMult.makeFootprint(ellipse)
    printFp(fp, "ellipseFp")

    del ellipse
    del core
    del point
    del fp

    core = geomEllipses.Axes(15, 15, 1.30)
    point = afwGeom.makePointD(0,0)
    ellipse = geomEllipses.AxesEllipse(core, point)
    fp = measMult.makeFootprint(ellipse)
    printFp(fp, "circleFp")

    del core
    del ellipse
    del point
    del fp

if __name__ == '__main__':
    main()
    

