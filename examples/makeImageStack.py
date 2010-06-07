import sys
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.afw.detection as afwDet
import lsst.meas.algorithms as measAlg
import lsst.meas.multifit as measMult
import numpy
import numpy.random
import lsst.afw.display.ds9 as ds9

def makeImageStack(model, depth, ra, dec):
    psf = measAlg.createPSF("DoubleGaussian", 7, 7, 1.0)
    crVal = afwGeom.makePointD(ra,dec)
    crPix = afwGeom.makePointD(0,0)    
    wcs = afwImage.createWcs(crVal, crPix, 1., 0., 0., 1.)

    fp = model.computeProjectionFootprint(psf, wcs)
    projection = model.makeProjection(psf, wcs, fp)

    bbox = fp.getBBox()
    nPix = fp.getNpix();

    #print >> sys.stderr, "fp bbox: %s"%bbox
    modelImage = afwImage.MaskedImageF(bbox.getWidth(), bbox.getHeight())
    modelImage.setXY0(bbox.getX0(), bbox.getY0())

    imageVector = projection.computeModelImage()
    varianceVector = numpy.zeros(nPix, dtype=numpy.float32)

    measMult.expandImageF(fp, modelImage, imageVector, varianceVector)
    bitMask = afwImage.MaskU.addMaskPlane("MULTIFIT")
    numpy.set_printoptions(threshold=100000)
    #pirint >> sys.stderr, "imageVector: %s"%imageVector
    #print >> sys.stderr, "dLinear: %s"%projection.computeLinearParameterDerivative()
    #print >> sys.stderr, "dNonlinear: %s"%projection.computeNonlinearParameterDerivative()
    sigma = 0.5
    sigmaSq = sigma**2

    expList = measMult.CharacterizedExposureListF()
    afwRandom = afwMath.Random()

    for i in xrange(depth):
        #create an exposure whose size matches footprint size
        exp = measMult.CharacterizedExposureF( 
            bbox.getWidth(), bbox.getHeight(), wcs, psf
        )
        mi = exp.getMaskedImage()        
        mi.setXY0(bbox.getX0(), bbox.getY0())

        img = mi.getImage()
        afwMath.randomUniformImage(img, afwRandom)

        img += modelImage.getImage()
        del img

        var = mi.getVariance()
        var.set(sigmaSq)
        del var

        msk = mi.getMask()        
        afwDet.setMaskFromFootprint(msk, fp, bitMask)
        msk^= bitMask
        del msk

        
        expList.append(exp)

    return expList

def writeImageStack(depth, baseName):
    flux = 1.0
    position = afwGeom.makePointD(45,45)
    psModel = measMult.createPointSourceModel(flux, position)
    
    filenameFrmt = basename + "_%%0%dd.fits"%len(str(depth-1))
    expList = makeImageStack(psModel, depth, position.getX(), position.getY())
    for i, exp in enumerate(expList):
        filename = filenameFrmt%i
        exp.writeFits(filename)

def displayImageStack(depth):
    flux = 1.0
    position = afwGeom.makePointD(45,45)
    psModel = measmultifit.createPointSourceModel(flux, position)
    expList = makeImageStack(psModel, depth, position.getX(), position.getY())

    for i, exp in enumerate(expList):
        ds9.mtv(exp, title=str(i))

if __name__ == "__main__":
    argc = len(sys.argv)
    if argc == 1:
        displayImageStack(1)
    else if argc == 2:
        displayImageStack(int(sys.argv[1]))
    else if argc == 3:
        writeImageStack(int(sys.argv[1]), sys.argv[2])
    else:      
        usage =
        """
        Usage:
            makeImageStack [<depth> [<output prefix>]]

        <depth> - number of images to produce [1, ) defaults to 1
        <output prefix> - specify the prefix of the output path.
                  If an output prefix is present, the image stack will be 
                  persisted to .fits files. Each file will be named:
                  "prefix_n.fits" where n is an integer from 0 to depth
        """

