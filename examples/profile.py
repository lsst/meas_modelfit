import lsst.meas.multifit
import lsst.afw.geom.ellipses
import lsst.afw.detection
import numpy
from matplotlib import pyplot

def bin(array, n=4):
    r = numpy.zeros(array.size / n, dtype=float)
    for i in range(n):
        r += array[i:array.size-n+i:n]
    return r / n
    

def main(basis):
    ellipse = lsst.afw.geom.ellipses.Ellipse(
        lsst.afw.geom.ellipses.Axes(12, 10, 0.25),
        lsst.afw.geom.Point2D(0.0, 0.0)
        )
    bounds = lsst.afw.geom.ellipses.Ellipse(ellipse)
    bounds.scale(5.0)
    footprint = lsst.afw.detection.Footprint(bounds)
    ibox = lsst.afw.geom.BoxI(footprint.getBBox())
    dbox = lsst.afw.geom.BoxD(ibox)
    extent = (dbox.getMinX(), dbox.getMaxX(), dbox.getMinY(), dbox.getMaxY())
    modelMatrix = numpy.zeros((footprint.getArea(), basis.getSize()), dtype=float)
    modelImages = numpy.zeros((ibox.getHeight(), ibox.getWidth(), basis.getSize()), dtype=float)
    basis.evaluate(modelMatrix, footprint, ellipse)
    lsst.afw.detection.expandArray(footprint, modelMatrix, modelImages, ibox.getMin())
    radii = numpy.linspace(1E-3, 10, 100)
    profileMatrix = numpy.zeros((radii.size, basis.getSize()), dtype=float)
    basis.evaluateRadialProfile(profileMatrix, radii)
    elements = [n for n in range(basis.getSize()) if (numpy.abs(profileMatrix[:,n]) > 1E-15).any()]
    gt = ellipse.getGridTransform()
    xg, yg = numpy.meshgrid(
        range(ibox.getMinX(), ibox.getMaxX() + 1),
        range(ibox.getMinY(), ibox.getMaxY() + 1)
        )
    xgt = gt[0, 0] * xg + gt[0, 1] * yg + gt[0, 2]
    ygt = gt[1, 0] * xg + gt[1, 1] * yg + gt[1, 2]
    rgt = (xgt**2 + ygt**2)**0.5
    for i, n in enumerate(elements):
        pyplot.subplot(2, len(elements), i + 1)
        pyplot.imshow(modelImages[:,:,n], origin='lower', interpolation='nearest', extent=extent)
        pyplot.title(n)
        pyplot.subplot(2, len(elements), len(elements) + i + 1)
        pyplot.plot(radii, profileMatrix[:,n])
        pyplot.plot(bin(rgt.ravel()), bin(modelImages[:,:,n].ravel()), ',', alpha=0.1)
    pyplot.show()
