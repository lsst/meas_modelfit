import lsst.meas.multifit
import lsst.afw.geom.ellipses
import lsst.afw.detection
import lsst.afw.math.shapelets
import numpy
from matplotlib import pyplot

import lsst.utils.tests as utilsTests
import unittest

def makeBasisImages(basis, ellipse, factor):
    n = basis.getSize()
    bounds = lsst.afw.geom.ellipses.Ellipse(ellipse)
    bounds.scale(factor)
    box = lsst.afw.geom.Box2I(bounds.computeEnvelope())
    footprint = lsst.afw.detection.Footprint(box)
    array = numpy.zeros((footprint.getArea(), basis.getSize()), dtype=float)
    basis.evaluate(array, footprint, ellipse)
    return array.reshape(box.getHeight(), box.getWidth(), n)

class IntegrationTestCase(unittest.TestCase):

    def compare(self, basis, ellipse, factor):
        images = makeBasisImages(basis, ellipse, factor)
        i1 = images.sum(axis=0).sum(axis=0)
        i2 = numpy.zeros(basis.getSize(), dtype=float)
        basis.integrate(i2)
        self.assert_(numpy.allclose(i1, i2, rtol=1E-2, atol=1E-1), "%s\n!=\n%s" % (i1, i2))

    def testIntegration(self):
        basis1 = lsst.meas.multifit.ShapeletModelBasis.make(2, 1.2)
        basis2 = lsst.meas.multifit.SourceMeasurement.loadBasis("exponential")
        psf = lsst.afw.math.shapelets.ShapeletFunction(
            2, lsst.afw.math.shapelets.HERMITE,
            lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes(1.2, 1.0, 0.0))
            )
        psf.getCoefficients()[:] = numpy.random.randn(psf.getCoefficients().size)
        psf.getCoefficients()[0] = 5.0
        psf.normalize()
        convolved1 = basis1.convolve(psf)
        convolved2 = basis2.convolve(psf)
        ellipse = lsst.afw.geom.ellipses.Ellipse(lsst.afw.geom.ellipses.Axes(8.0, 7.0, 0.0))
        self.compare(basis1, ellipse, 5)
        self.compare(basis2, ellipse, 5)
        self.compare(convolved1, ellipse, 5)
        self.compare(convolved2, ellipse, 5)

def suite():
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(IntegrationTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)


def run(exit=False):
    """Run the tests"""
    utilsTests.run(suite(), exit)

if __name__ == "__main__":
    run(True)
