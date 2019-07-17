#
# LSST Data Management System
#
# Copyright 2008-2016  AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <https://www.lsstcorp.org/LegalNotices/>.
#
import os
import unittest
import numpy

import lsst.utils.tests
import lsst.meas.modelfit

try:
    import scipy.stats
except ImportError:
    scipy = None


class MixtureTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        numpy.random.seed(500)
        self.rng = lsst.afw.math.Random("MT19937", 500)

    @staticmethod
    def makeRandomMixture(nDim, nComponents, df=float("inf")):
        componentList = []
        for i in range(nComponents):
            mu = numpy.random.randn(nDim)*4
            a = numpy.random.randn(nDim+1, nDim)
            sigma = numpy.dot(a.transpose(), a) + numpy.identity(nDim)
            componentList.append(lsst.meas.modelfit.Mixture.Component(numpy.random.rand(), mu, sigma))
        return lsst.meas.modelfit.Mixture(nDim, componentList, df)

    def testWrappers(self):
        """Test that we correctly wrapped tricky things.
        """
        l1 = []
        l1.append(lsst.meas.modelfit.MixtureComponent(1))
        l1.append(lsst.meas.modelfit.MixtureComponent(1))
        l1.append(lsst.meas.modelfit.MixtureComponent(1))
        l1[0].weight = 1.0
        l1[0].setMu(numpy.array([1.0], dtype=float))
        l1[0].setSigma(numpy.array([[4.0]], dtype=float))
        l1[1].weight = 0.5
        l1[2].weight = 0.5
        m1 = lsst.meas.modelfit.Mixture(1, l1)
        self.assertEqual(m1[0].weight, 0.5)
        self.assertEqual([0.5, 0.25, 0.25], [c.weight for c in m1])
        self.assertFloatsAlmostEqual(m1[0].getMu(), numpy.array([1.0], dtype=float))
        self.assertFloatsAlmostEqual(m1[0].getSigma(), numpy.array([4.0], dtype=float))
        self.assertFloatsAlmostEqual(m1.evaluate(m1[1], numpy.array([0.0], dtype=float)),
                                     m1[1].weight*(2.0*numpy.pi)**(-0.5))
        self.assertFloatsAlmostEqual(m1.evaluate(numpy.array([0.0], dtype=float)),
                                     (m1[0].weight*numpy.exp(-0.125)/2 + m1[1].weight + m1[2].weight) *
                                     (2.0*numpy.pi)**(-0.5))

    def testGaussian(self):
        """Test that our implementations for a single-component Gaussian are correct.
        """
        m = self.makeRandomMixture(2, 1)
        mu = m[0].getMu()
        sigma = m[0].getSigma()
        fisher = numpy.linalg.inv(sigma)
        x = numpy.random.randn(20, 2)
        p = numpy.zeros(20, dtype=float)
        m.evaluate(x, p)
        z = ((x - mu)[:, numpy.newaxis, :] * fisher[numpy.newaxis, :, :, ] *
             (x - mu)[:, :, numpy.newaxis]).sum(axis=2).sum(axis=1)
        self.assertFloatsAlmostEqual(p, numpy.exp(-0.5*z) / numpy.linalg.det(2*numpy.pi*sigma)**0.5)
        x = numpy.zeros((1000000, 2), dtype=float)
        m.draw(self.rng, x)
        self.assertFloatsAlmostEqual(x.mean(axis=0), mu, rtol=2E-2)
        self.assertFloatsAlmostEqual(numpy.cov(x, rowvar=False), sigma, rtol=3E-2)

    @unittest.skipIf(scipy is None, "Test requires SciPy")
    def testGaussianSciPy(self):
        m = self.makeRandomMixture(2, 1)
        x = numpy.zeros((1000000, 2), dtype=float)
        m.draw(self.rng, x)
        self.assertGreater(scipy.stats.normaltest(x[:, 0])[1], 0.05)
        self.assertGreater(scipy.stats.normaltest(x[:, 1])[1], 0.05)

    @unittest.skipIf(scipy is None, "Test requires SciPy")
    def testStudentsT(self):
        """Test that our implementations for a single-component Student's T are correct.
        """
        for df in [4, 8]:
            m = self.makeRandomMixture(1, 1, df=df)
            mu = m[0].getMu()
            sigma = m[0].getSigma()
            x = numpy.random.randn(20, 1)
            p = numpy.zeros(20, dtype=float)
            m.evaluate(x, p)
            x = x.reshape(20)
            z = (x - mu)/(sigma**0.5)
            self.assertFloatsAlmostEqual(p, scipy.stats.t.pdf(z, df)/sigma**0.5)
            x = numpy.zeros((1000000, 1), dtype=float)
            m.draw(self.rng, x)
            self.assertFloatsAlmostEqual(x.mean(), mu, rtol=5E-2)
            self.assertFloatsAlmostEqual(x.var(), sigma * df / (df - 2), rtol=5E-2)
            self.assertLess(scipy.stats.normaltest(x)[1], 0.05)

    def testPersistence(self):
        """Test table-based persistence of Mixtures"""
        filename = "testMixturePersistence.fits"
        mix1 = self.makeRandomMixture(3, 4, df=3.5)
        mix1.writeFits(filename)
        mix2 = lsst.meas.modelfit.Mixture.readFits(filename)
        self.assertEqual(mix1.getDegreesOfFreedom(), mix2.getDegreesOfFreedom())
        self.assertEqual(len(mix1), len(mix2))
        for c1, c2 in zip(mix1, mix2):
            self.assertFloatsAlmostEqual(c1.weight, c2.weight)
            self.assertFloatsAlmostEqual(c1.getMu(), c2.getMu())
            self.assertFloatsAlmostEqual(c1.getSigma(), c2.getSigma())
        os.remove(filename)

    def testDerivatives(self):
        epsilon = 1E-7
        g = self.makeRandomMixture(3, 4)
        t = self.makeRandomMixture(4, 3, df=4.0)

        def doTest(mixture, point):
            n = mixture.getDimension()
            # Compute numeric first derivatives
            testPoints = numpy.zeros((2*n, n), dtype=float)
            testPoints[:, :] = point[numpy.newaxis, :]
            for i in range(n):
                testPoints[i, i] += epsilon
                testPoints[n+i, i] -= epsilon
            testValues = numpy.zeros(2*n, dtype=float)
            mixture.evaluate(testPoints, testValues)
            numericGradient = numpy.zeros(n, dtype=float)
            for i in range(n):
                numericGradient[i] = (testValues[i] - testValues[n+i]) / (2.0 * epsilon)
            # Compute numeric second derivatives from analytic first derivatives
            numericHessian = numpy.zeros((n, n), dtype=float)
            testGrad1 = numpy.zeros(n, dtype=float)
            testGrad2 = numpy.zeros(n, dtype=float)
            testHessian = numpy.zeros((n, n), dtype=float)
            for i in range(n):
                testPoint = point.copy()
                testPoint[i] += epsilon
                mixture.evaluateDerivatives(testPoint, testGrad1, testHessian)
                testPoint[i] -= 2.0*epsilon
                mixture.evaluateDerivatives(testPoint, testGrad2, testHessian)
                numericHessian[i, :] = (testGrad1 - testGrad2) / (2.0 * epsilon)
            # Compute analytic derivatives and compare
            analyticGradient = numpy.zeros(n, dtype=float)
            analyticHessian = numpy.zeros((n, n), dtype=float)
            mixture.evaluateDerivatives(point, analyticGradient, analyticHessian)
            self.assertFloatsAlmostEqual(analyticGradient, numericGradient, rtol=1.5E-6)
            self.assertFloatsAlmostEqual(analyticHessian, numericHessian, rtol=1E-6)

        for x in numpy.random.randn(10, g.getDimension()):
            doTest(g, x)

        for x in numpy.random.randn(10, t.getDimension()):
            doTest(t, x)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
