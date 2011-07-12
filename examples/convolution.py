import lsst.meas.multifit
import lsst.afw.geom.ellipses
import lsst.afw.detection
import lsst.afw.math.shapelets
import numpy
from matplotlib import pyplot

def makeBasisImages(basis, ellipse, factor):
    n = basis.getSize()
    bounds = lsst.afw.geom.ellipses.Ellipse(ellipse)
    bounds.scale(factor)
    box = lsst.afw.geom.Box2I(bounds.computeEnvelope())
    footprint = lsst.afw.detection.Footprint(box)
    array = numpy.zeros((footprint.getArea(), basis.getSize()), dtype=float)
    basis.evaluate(array, footprint, ellipse)
    return array.reshape(box.getHeight(), box.getWidth(), n)

def plotBasisImages(images):
    pyplot.figure()
    for n in range(images.shape[2]):
        pyplot.subplot(1, images.shape[2], n + 1)
        pyplot.imshow(images[:,:,n])
        pyplot.axis("off")
        pyplot.title("%0.04f" % images[:,:,n].sum(), fontsize=10)

def main():
    numpy.set_printoptions(suppress=True, linewidth=180)
    #basis = lsst.meas.multifit.ShapeletModelBasis.make(2, 1.2)
    basis = lsst.meas.multifit.SourceMeasurement.loadBasis("exponential")
    psf = lsst.afw.math.shapelets.ShapeletFunction(
        2, lsst.afw.math.shapelets.HERMITE,
        lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes(2.0, 1.8, 0.0))
        )
    psf.getCoefficients()[:] = numpy.random.randn(psf.getCoefficients().size)
    psf.getCoefficients()[0] = 5.0
    psf.normalize()
    convolved = basis.convolve(psf)
    ellipse1 = lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes(8.0, 8.0, 0.0))
    ellipse2 = lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes(12.0, 12.0, 0.0))
    images1 = makeBasisImages(basis, ellipse1, 10)
    images2 = makeBasisImages(basis, ellipse2, 10)
    images1c = makeBasisImages(convolved, ellipse1, 10)
    images2c = makeBasisImages(convolved, ellipse2, 10)
    plotBasisImages(images1)
    plotBasisImages(images2)
    plotBasisImages(images1c)
    plotBasisImages(images2c)
    pyplot.show()
    
if __name__ == "__main__":
    main()
