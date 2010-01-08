import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.meas.algorithms as measAlg
import lsst.meas.multifit as measMult
import numpy

def makeImageStack(depth):
    psFactory = measMult.PointSourceModelFactory()
    psModel = psFactory.makeModel(1.0, afwGeom.PointD(0,0))

    
    psf = measAlg.createPSF("DoubleGaussian", 9, 9, 3)
    wcs = afwImage.createWcs(
        afwImage.PointD(0,0), 
        afwImage.POintD(0,0), 
        1, 0, 0, 1)
    
    fp = psModel.computeProjectionFootprint(psf, wcs)
    psProjection = psModel.makeProjection(psf, wcs, fp)

    bbox = fp.getBBox()
    
    exp = makeExposureD(bbox.getWidth(), bbox.getHeight())
    mi = exp.getMaskedImage()
    mi.setXY0(bbox.getX0(), bbox.getY0())

    epxList []

        
    
