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

def compare(basis, ellipse, factor):
    images = makeBasisImages(basis, ellipse, factor)
    i1 = images.sum(axis=0).sum(axis=0)
    i2 = numpy.zeros(basis.getSize(), dtype=float)
    basis.integrate(i2)
    i2 *= ellipse.getCore().getArea() / numpy.pi
    print i1
    print i2

def main():
    numpy.set_printoptions(suppress=True, linewidth=180)
    basis1 = lsst.meas.multifit.ShapeletModelBasis.make(2, 1.2)
    basis2 = lsst.meas.multifit.SourceMeasurement.loadBasis(8)
    psf = lsst.afw.math.shapelets.ShapeletFunction(
        2, lsst.afw.math.shapelets.HERMITE,
        lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes(2.0, 1.8, 0.0))
        )
    psf.getCoefficients()[:] = numpy.random.randn(psf.getCoefficients().size)
    psf.getCoefficients()[0] = 5.0
    psf.normalize()
    convolved1 = basis1.convolve(psf)
    convolved2 = basis2.convolve(psf)
    ellipse = lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes(5.0, 5.0, 0.0))
    compare(basis1, ellipse, 10)
    compare(basis2, ellipse, 20)
    
if __name__ == "__main__":
    main()
