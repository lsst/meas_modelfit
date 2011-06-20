import lsst.afw.math.shapelets
import lsst.afw.detection
import lsst.afw.geom.ellipses
import lsst.meas.multifit.utils
import numpy
try:
    from matplotlib import pyplot
except ImportError:
    pass

def compare(basis, ellipse, factor=15, display=False):
    n = basis.getSize()
    bounds = lsst.afw.geom.ellipses.Ellipse(ellipse)
    bounds.scale(factor)
    box = lsst.afw.geom.Box2I(bounds.computeEnvelope())
    footprint = lsst.afw.detection.Footprint(box)
    array = numpy.zeros((footprint.getArea(), basis.getSize()), dtype=float)
    basis.evaluate(array, footprint, ellipse)
    shape = (box.getHeight(), box.getWidth())

    m1 = basis.computeInnerProductMatrix(ellipse.getCore())
    m2 = numpy.zeros_like(m1)

    if display:
        pyplot.figure(figsize=(12, 2))
        for i in range(n):
            axes = pyplot.subplot(1, n, i + 1)
            axes.imshow(array[:,i].reshape(*shape))
            axes.axis("off")
        pyplot.figure(figsize=(10, 10))
    for i in range(n):
        for j in range(i + 1):
            image = (array[:,i] * array[:,j]).reshape(*shape)
            m2[i,j] = image.sum()
            m2[j,i] = m2[i,j]
            if display:
                axes = pyplot.subplot(n, n, i*n + j + 1)
                axes.imshow(image)
                axes.axis("off")
                axes = pyplot.subplot(n, n, j*n + i + 1)
                axes.imshow(image)
                axes.axis("off")

    print m1
    print m2
    mask = numpy.abs(m2) > 1E-8
    print m2[mask] / m1[mask]

def main(display=False):
    numpy.set_printoptions(suppress=True, linewidth=180)
    basis1 = lsst.meas.multifit.ShapeletModelBasis.make(2, 1.2)
    basis2 = lsst.meas.multifit.utils.loadBasis("ed+06:2000")
    psf = lsst.afw.math.shapelets.ShapeletFunction(
        2, lsst.afw.math.shapelets.HERMITE,
        lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes(2.0, 1.8, 0.0))
        )
    psf.getCoefficients()[:] = numpy.random.randn(psf.getCoefficients().size)
    convolved1 = basis1.convolve(psf)
    convolved2 = basis2.convolve(psf)
    ellipse = lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes(5.0, 4.0, 0.0))
    #compare(basis1, ellipse, display=display)
    compare(convolved1, ellipse, display=display)
    #compare(basis2, ellipse, display=display)
    compare(convolved2, ellipse, display=display)
    if display:
        pyplot.show()

if __name__ == "__main__":
    main(False)

