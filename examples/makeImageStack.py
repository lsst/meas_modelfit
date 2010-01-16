import sys
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDet
import lsst.meas.algorithms as measAlg
import lsst.meas.multifit as measMult
import numpy
import numpy.random

def makeImageStack(model, depth, ra, dec):
    psf = measAlg.createPSF("DoubleGaussian", 19, 19, 2)
    crVal = afwImage.PointD(ra,dec)
    crPix = afwImage.PointD(0,0)    
    wcs = afwImage.createWcs(crVal, crPix, 0.0001, 0, 0, 0.0001)

    fp = model.computeProjectionFootprint(psf, wcs)
    projection = model.makeProjection(psf, wcs, fp)

    bbox = fp.getBBox()
    nPix = fp.getNpix();

    print "fp bbox: %s"%bbox
    imageVector = projection.computeModelImage()

    numpy.set_printoptions(threshold=100000)
    print "imageVector: %s"%imageVector
    print "dLinear: %s"%projection.computeLinearParameterDerivative()
    print "dNonlinear: %s"%projection.computeNonlinearParameterDerivative()
    sigma = 0.5

    expList = measMult.CharacterizedExposureListD()
    for i in xrange(depth):
        #create an exposure whose size matches footprint size
        exp = measMult.CharacterizedExposureD( 
            bbox.getWidth(), bbox.getHeight(), wcs, psf
        )
        mi = exp.getMaskedImage()
        mi.setXY0(bbox.getX0(), bbox.getY0())
        
        #set image to be model projection's modelImage
        #set variance to be a flat field
        varianceVector = numpy.zeros(nPix, dtype=numpy.float32)
        varianceVector[:] = sigma**2
        randomNormal = numpy.random.normal(scale=sigma, size=nPix)
        noisyImage = imageVector + randomNormal
        measMult.expandImageD(fp, mi, noisyImage, varianceVector)

        #set the mask to be the inverse of the footprint
        mask = mi.getMask()
        mask = 1
        afwDet.setMaskFromFootprint(mi.getMask(), fp, 0)

        expList.append(exp)

    return expList

def writeImageStack(depth, baseName):
    psFactory = measMult.PointSourceModelFactory()
    psModel = psFactory.makeModel(1.0, afwGeom.makePointD(45,45))
    
    expList = makeImageStack(psModel, depth, 45, 45)
    for i, exp in enumerate(expList):
        filename = baseName +"_%d"%i
        exp.writeFits(filename)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Two arguments required: number of images, and base name for output files"
    writeImageStack(int(sys.argv[1]), sys.argv[2])





    
