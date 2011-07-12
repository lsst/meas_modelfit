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
    xr = numpy.arange(box.getBeginX(), box.getEndX(), 1)
    yr = numpy.arange(box.getBeginY(), box.getEndY(), 1)
    images = array.reshape(box.getHeight(), box.getWidth(), n)
    xg, yg = numpy.meshgrid(xr, yr)
    return images, xg, yg

def transformCoords(ellipse, x, y):
    t = ellipse.getGridTransform().getMatrix()
    xt = x * t[0,0] + y * t[0,1] + t[0,2]
    yt = x * t[1,0] + y * t[1,1] + t[1,2]
    return xt, yt

def measureImageMoments(images, x, y):
    if len(images.shape) == 2:
        images = images.reshape(images.shape[0], images.shape[1], 1)
    multipoles = numpy.zeros((6, images.shape[2]), dtype=float)
    for i in range(images.shape[2]):
        multipoles[0,i] = images[:,:,i].sum()
        multipoles[1,i] = (x * images[:,:,i]).sum()
        multipoles[2,i] = (y * images[:,:,i]).sum()
        multipoles[3,i] = (x**2 * images[:,:,i]).sum()
        multipoles[4,i] = (y**2 * images[:,:,i]).sum()
        multipoles[5,i] = (x * y * images[:,:,i]).sum()
    return multipoles

def makeMomentsEllipse(m):
    m = m.copy()
    m[1:] /= m[0]
    m[3] -= m[1]**2
    m[4] -= m[2]**2
    m[5] -= m[1] * m[2]
    return lsst.afw.geom.ellipses.Ellipse(
        lsst.afw.geom.ellipses.Quadrupole(m[3], m[4], m[5]),
        lsst.afw.geom.Point2D(m[1], m[2])
        )

def plotBasisImages(images):
    pyplot.figure()
    for n in range(images.shape[2]):
        pyplot.subplot(1, images.shape[2], n + 1)
        pyplot.imshow(images[:,:,n])
        pyplot.axis("off")

def main():
    numpy.set_printoptions(suppress=True, linewidth=180)
    basis = lsst.meas.multifit.ShapeletModelBasis.make(2, 1.5)
    #basis = lsst.meas.multifit.SourceMeasurement.loadBasis("exponential")
    ellipse1 = lsst.afw.geom.ellipses.Ellipse(
        lsst.afw.geom.ellipses.Quadrupole(lsst.afw.geom.ellipses.Axes(8.0, 6.0, 0.0))
        )
    ellipse2 = lsst.afw.geom.ellipses.Ellipse(
        lsst.afw.geom.ellipses.Quadrupole(lsst.afw.geom.ellipses.Axes(12.0, 10.0, 0.45))
        )
    images1, x1, y1 = makeBasisImages(basis, ellipse1, 10)
    images2, x2, y2 = makeBasisImages(basis, ellipse2, 10)
    xt1, yt1 = transformCoords(ellipse1, x1, y1)
    xt2, yt2 = transformCoords(ellipse2, x2, y2)
    multipoles1a = measureImageMoments(images1, xt1, yt1)
    multipoles2a = measureImageMoments(images2, xt2, yt2)
    mm = basis.getMultipoleMatrix()
    mm1b = lsst.meas.multifit.MultipoleMatrix(measureImageMoments(images1, x1, y1))
    mm2b = lsst.meas.multifit.MultipoleMatrix(measureImageMoments(images2, x2, y2))
    coeff1 = numpy.random.randn(basis.getSize())
    coeff2 = numpy.random.randn(basis.getSize())
    m1b = numpy.dot(mm1b.getArray(), coeff1)
    m2b = numpy.dot(mm2b.getArray(), coeff2)
    mm1b.transform(ellipse1.getGridTransform())
    mm2b.transform(ellipse2.getGridTransform())
    #plotBasisImages(images1)
    #plotBasisImages(images2)
    print numpy.allclose(multipoles1a, multipoles2a, rtol=1E-4, atol=1E-4)
    print numpy.allclose(multipoles1a, mm.getArray(), rtol=1E-4, atol=1E-4)
    print numpy.allclose(multipoles1a, mm1b.getArray(), rtol=1E-4, atol=1E-4)
    print numpy.allclose(multipoles1a, mm2b.getArray(), rtol=1E-4, atol=1E-4)
    ellipse1a = lsst.afw.geom.ellipses.Ellipse(ellipse1)
    ellipse2a = lsst.afw.geom.ellipses.Ellipse(ellipse2)
    mm.applyMoments(ellipse1a, coeff1)
    mm.applyMoments(ellipse2a, coeff2)
    print ellipse1a
    print makeMomentsEllipse(m1b)
    print ellipse2a
    print makeMomentsEllipse(m2b)
    #pyplot.show()
    
if __name__ == "__main__":
    main()
