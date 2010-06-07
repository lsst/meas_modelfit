import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImg
import lsst.afw.detection as afwDet
import lsst.meas.multifit as measMult
import numpy


def main():
    fp = afwDet.Footprint(afwImg.BBox(afwImg.PointI(0,0), 2, 3))
    fp.normalize()
    window = afwGeom.BoxI(afwGeom.makePointI(1,0), afwGeom.makeExtentI(1,3))
    wfp = measMult.WindowedFootprint(fp, window)
    full = numpy.array([[1., 2.], [3., 4.], [5., 6.]])
    compressed = numpy.zeros(fp.getNpix())
    wfp.compress(full, compressed)

    print full
    print window
    print compressed

if __name__ == '__main__':
    main()
    

