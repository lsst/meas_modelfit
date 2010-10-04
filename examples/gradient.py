import numpy
import lsst.meas.multifit as measMult
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as geomEllipses
import lsst.afw.math as afwMath
import lsst.afw.detection as afwDet
import math

def makeModelImage(model, psf, affine, noiseFactor):
    fp = model.computeProjectionFootprint(psf, affine)
    nPix = fp.getNpix()
    box = fp.getBBox()
    proj = model.makeProjection(psf, affine, fp)
    modelImage = afwImage.MaskedImageF(box.getWidth(), box.getHeight())
    modelImage.setXY0(box.getX0(), box.getY0())
    
    imageVector = proj.computeModelImage()
    varianceVector = numpy.zeros_like(imageVector)

    measMult.expandImageF(fp, modelImage, imageVector, varianceVector)

    afwRandom = afwMath.Random()
    randomImg = afwImage.ImageF(modelImage.getDimensions())
    afwMath.randomGaussianImage(randomImg, afwRandom)
    randomImg *= noiseFactor
    img = modelImage.getImage()
    img += randomImg

    stats = afwMath.makeStatistics(modelImage, afwMath.VARIANCE)
    variance = stats.getValue(afwMath.VARIANCE)
    #variance = 0.25
    modelImage.getVariance().set(variance)
    return modelImage

def makeModelExposure(model, psf, wcs, noiseFactor=0):
    affine = wcs.linearizeSkyToPixel(model.getPosition())
    modelImage = makeModelImage(model, psf, affine, noiseFactor)
    
    exp = afwImage.ExposureF(modelImage, wcs)
    exp.setPsf(psf)
    return exp

def computeValue(modelEvaluator):
    modeled = modelEvaluator.computeModelImage()
    measured = modelEvaluator.getWeightedData()
    residual = measured-modeled

    return 0.5*residual.T*residual

def computeGradient(modelEvaluator, checkNumericGradient):
    modeled = modelEvaluator.computeModelImage()
    measured = modelEvaluator.getWeightedData()

    residual = measured - modeled
    lpd = modelEvaluator.computeLinearParameterDerivative()
    npd = modelEvaluator.computeNonlinearParameterDerivative()
   
    gradLinear = (-lpd.T)*residual
    gradNonlinear = (-npd.T)*residual

    grad = numpy.concatenate((gradLinear, gradNonlinear))

    if(checkNumericGradient):
        nonlinear = modelEvaluator.getNonlinearParameters()
        linear = modelEvaluator.getLinearParameters()
    
        h = [1e-4, 1e-4, 1e-4]
        
        columnVectors = []
        i = 0
        for n in range(len(linear)):
            linear[n] += h[i] 
            
            modelEvaluator.setLinearParameters(linear)
            plus = modelEvaluator.computeModelImage()
            
            linear[n] -= 2*h[i]
            modelEvaluator.setLinearParameters(linear)
            minus = modelEvaluator.computeModelImage()

            partial = (plus - minus) / (2*h[i])
            columnVectors.append(partial)

            linear[n] += h[i]
            i+=1

        for n in range(len(nonlinear)):
            nonlinear[n] += h[i]        
            modelEvaluator.setNonlinearParameters(nonlinear)
            plus = modelEvaluator.computeModelImage()

            nonlinear[n] -= 2*h[i]

            modelEvaluator.setNonlinearParameters(nonlinear)
            minus = modelEvaluator.computeModelImage()
            partial = (plus - minus)/(2*h[i])
        
            columnVectors.append(partial)
            nonlinear[n] += h[i]
            i+=1
        numeric = numpy.concatenate(columnVectors, 1)

        print numeric.shape
        numericGradient = (-numeric.T)*residual
        print "numeric gradient", numericGradient

    return grad


def main():
    numpy.set_printoptions(threshold=numpy.nan)

    centroid = afwGeom.makePointD(45,45)
    psf = afwDet.createPsf("DoubleGaussian", 9, 9, 1.0)
    wcs = afwImage.createWcs(centroid, afwGeom.makePointD(0,0), 0.0001, 0., 0., 0.0001)
    #affine = wcs.linearizePixelToSky(centroid)
    affine = afwGeom.AffineTransform()

    flux = 20.0
    axes = geomEllipses.Axes(25,30,0).transform(affine)
    sersicIndex = 1.5

    model = measMult.createPointSourceModel(flux, centroid)
    #model = measMult.createExponentialModel(flux, centroid, axes)

    modelImage = makeModelImage(model, psf, affine, 4.0)
    eval = measMult.ModelEvaluator(model, affine)
    exp = makeModelExposure(model, psf, wcs)
    expList = [exp]

    images = measMult.MaskedImageListF()
    images.append(modelImage)
    psfs = measMult.PsfList()
    psfs.append(psf)
    transforms = measMult.TransformList()
    transforms.append(affine)

    eval.setData( images, psfs, transforms)
    #eval.setExposures(expList)

    print "linear:", eval.getLinearParameters()
    print "nonlinear:", eval.getNonlinearParameters()
    print computeValue(eval)
    print computeGradient(eval, True)

if __name__ == '__main__':
    main()
