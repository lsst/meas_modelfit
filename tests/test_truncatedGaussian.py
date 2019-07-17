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
import unittest
import numpy
try:
    import scipy.integrate
    import scipy.stats
    import scipy.special
except ImportError:
    scipy = None

import lsst.log
import lsst.log.utils
import lsst.utils.tests
import lsst.meas.modelfit


if False:
    lsst.log.utils.traceSetAt("meas.modelfit.integrals", 5)
    lsst.log.utils.traceSetAt("meas.modelfit.TruncatedGaussian", 5)


class TruncatedGaussianTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        numpy.random.seed(500)

    def check1d(self, mu, hessian, tg):
        evaluator = tg.evaluate()
        logEvaluator = tg.evaluateLog()
        dist = scipy.stats.norm(loc=mu[0], scale=hessian[0, 0]**-0.5)
        self.assertFloatsAlmostEqual(1.0 - dist.cdf(0.0), tg.getUntruncatedFraction())
        eps = 1E-7
        if numpy.all(mu >= 0.0):
            self.assertFloatsAlmostEqual(logEvaluator(mu), tg.getLogPeakAmplitude())
            self.assertGreater(logEvaluator(mu+eps), tg.getLogPeakAmplitude())
            self.assertGreater(logEvaluator(mu-eps), tg.getLogPeakAmplitude())
        peak = tg.maximize()
        self.assertGreater(evaluator(peak), 0.0)
        self.assertLess(evaluator(peak+eps), evaluator(peak))
        self.assertLess(evaluator(peak-eps), evaluator(peak))

        def altLogEval(x):
            if numpy.any(x < 0):
                return float("inf")
            return tg.getLogPeakAmplitude() + 0.5*hessian[0, 0]*(x-mu[0])**2
        for alpha in (numpy.random.randn(10, 1) * hessian[0, 0]**-0.5 + mu[0]):
            x1 = logEvaluator(alpha)
            x2 = altLogEval(alpha[0])
            if numpy.isfinite(x1) and numpy.isfinite(x2):
                self.assertFloatsAlmostEqual(x1, x2, rtol=1E-14)
            else:
                self.assertEqual(x1, x2)
        integral, check = self.integrate1d(tg)
        self.assertLess(check, 1E-7)
        self.assertFloatsAlmostEqual(integral, numpy.exp(-tg.getLogIntegral()), atol=check)

    def check2d(self, mu, hessian, tg, isDegenerate=False):
        evaluator = tg.evaluate()
        logEvaluator = tg.evaluateLog()
        unit1 = numpy.array([1.0, 0.0])
        unit2 = numpy.array([0.0, 1.0])
        eps = 1E-7
        if numpy.all(mu >= 0.0):
            self.assertFloatsAlmostEqual(logEvaluator(mu), tg.getLogPeakAmplitude())
            self.assertGreater(logEvaluator(mu + unit1*eps), tg.getLogPeakAmplitude())
            self.assertGreater(logEvaluator(mu - unit1*eps), tg.getLogPeakAmplitude())
            self.assertGreater(logEvaluator(mu + unit2*eps), tg.getLogPeakAmplitude())
            self.assertGreater(logEvaluator(mu - unit2*eps), tg.getLogPeakAmplitude())
        peak = tg.maximize()
        self.assertGreater(evaluator(peak), 0.0)
        self.assertLess(evaluator(peak + unit1*eps) / evaluator(peak), 1.0)
        self.assertLess(evaluator(peak - unit1*eps) / evaluator(peak), 1.0)
        self.assertLess(evaluator(peak + unit2*eps) / evaluator(peak), 1.0)
        self.assertLess(evaluator(peak - unit2*eps) / evaluator(peak), 1.0)

        def altLogEval(a):
            if numpy.any(a < 0):
                return float("inf")
            return tg.getLogPeakAmplitude() + 0.5*numpy.dot(numpy.dot(hessian, a - mu).transpose(), a - mu)
        for alpha in (numpy.random.randn(10, 2) * hessian.diagonal()**-0.5 + mu):
            x1 = logEvaluator(alpha)
            x2 = altLogEval(alpha)
            if numpy.isfinite(x1) and numpy.isfinite(x2):
                self.assertFloatsAlmostEqual(x1, x2, rtol=1E-14)
            else:
                self.assertEqual(x1, x2)
        integral, check = self.integrate2d(tg)
        self.assertLess(check, 1E-7)
        self.assertFloatsAlmostEqual(integral, numpy.exp(-tg.getLogIntegral()), atol=check)

    def integrate1d(self, tg):
        evaluator = tg.evaluate()

        def func(x):
            return evaluator(numpy.array([x]))
        return scipy.integrate.quad(func, 0.0, numpy.Inf)

    def integrate2d(self, tg):
        evaluator = tg.evaluate()

        def func(x, y):
            return evaluator(numpy.array([x, y]))
        return scipy.integrate.dblquad(func, 0.0, numpy.Inf, lambda x: 0.0, lambda x: numpy.Inf)

    @unittest.skipIf(scipy is None, "Test requires SciPy")
    def test1d(self):
        for i in range(5):
            sigma = (numpy.random.randn(1, 1)**2 + 1)*5
            mu = (numpy.random.randn(1))*3
            q0 = float(numpy.random.randn())
            hessian = numpy.linalg.inv(sigma)
            gradient = -numpy.dot(hessian, mu)
            tg1 = lsst.meas.modelfit.TruncatedGaussian.fromStandardParameters(mu, sigma)
            tg2 = lsst.meas.modelfit.TruncatedGaussian.fromSeriesParameters(q0, gradient, hessian)
            self.assertEqual(tg1.getLogIntegral(), 0.0)
            self.assertFloatsAlmostEqual(tg1.getLogPeakAmplitude(),
                                         (0.5*numpy.log(numpy.linalg.det(2.0*numpy.pi*sigma)) +
                                          numpy.log(tg1.getUntruncatedFraction())),
                                         rtol=1E-13)
            self.assertFloatsAlmostEqual(tg2.getLogPeakAmplitude(),
                                         q0 + 0.5*numpy.dot(mu, gradient),
                                         rtol=1E-13)
            self.check1d(mu, hessian, tg1)
            self.check1d(mu, hessian, tg2)

    @unittest.skipIf(scipy is None, "Test requires SciPy")
    def test2d(self):
        for i in range(5):
            x = numpy.linspace(-1, 1, 5)
            model = numpy.zeros((x.size, 2), dtype=float)
            model[:, 0] = x
            model[:, 1] = x**2 + x
            data = numpy.random.randn(x.size) + model[:, 0]*0.9 + model[:, 1]*1.1
            q0 = 0.5*float(numpy.dot(data, data))
            gradient = -numpy.dot(model.transpose(), data)
            hessian = numpy.dot(model.transpose(), model)
            sigma = numpy.linalg.inv(hessian)
            self.assertFloatsAlmostEqual(numpy.linalg.inv(sigma), hessian, rtol=1E-15, atol=1E-15)
            mu = -numpy.dot(sigma, gradient)
            tg1 = lsst.meas.modelfit.TruncatedGaussian.fromStandardParameters(mu, sigma)
            self.assertFloatsAlmostEqual(tg1.getLogPeakAmplitude(),
                                         (numpy.log(tg1.getUntruncatedFraction()) +
                                          0.5*numpy.log(numpy.linalg.det(2.0*numpy.pi*sigma))),
                                         rtol=1E-13)
            self.assertEqual(tg1.getLogIntegral(), 0.0)
            self.check2d(mu, hessian, tg1)
            tg2 = lsst.meas.modelfit.TruncatedGaussian.fromSeriesParameters(q0, gradient, hessian)
            self.assertFloatsAlmostEqual(tg2.getLogPeakAmplitude(),
                                         q0+0.5*numpy.dot(mu, gradient),
                                         rtol=1E-13)
            self.check2d(mu, hessian, tg2)

    @unittest.skipIf(scipy is None, "Test requires SciPy")
    def testDegenerate(self):
        for i in range(5):
            x = numpy.linspace(-1, 1, 5)
            model = numpy.zeros((x.size, 2), dtype=float)
            model[:, 0] = x
            model[:, 1] = 2*x
            data = numpy.random.randn(x.size) + model[:, 0]*0.9 + model[:, 1]*1.1
            q0 = 0.5*float(numpy.dot(data, data))
            gradient = -numpy.dot(model.transpose(), data)
            hessian = numpy.dot(model.transpose(), model)
            mu, _, _, _ = numpy.linalg.lstsq(model, data, rcond=None)
            tg = lsst.meas.modelfit.TruncatedGaussian.fromSeriesParameters(q0, gradient, hessian)
            self.assertFloatsAlmostEqual(tg.getLogPeakAmplitude(),
                                         q0+0.5*numpy.dot(mu, gradient),
                                         rtol=1E-13)
            self.check2d(mu, hessian, tg, isDegenerate=True)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
