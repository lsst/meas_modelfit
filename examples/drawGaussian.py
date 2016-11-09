#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2013 LSST Corporation.
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
# see <http://www.lsstcorp.org/LegalNotices/>.
#

from __future__ import print_function
from builtins import object
import numpy
import matplotlib

import lsst.pex.logging
import lsst.afw.math
import lsst.meas.modelfit.display

rng = lsst.afw.math.Random()
log = lsst.pex.logging.Debug("meas.modelfit.TruncatedGaussian", 10)


class TruncatedGaussianData(object):

    def __init__(self, tg, nSamples, strategy=None):
        self.tg = tg
        self.dimensions = ['x', 'y']
        self.values = numpy.zeros((nSamples, self.tg.getDim()), dtype=lsst.meas.modelfit.Scalar)
        self.weights = numpy.zeros(nSamples, dtype=lsst.meas.modelfit.Scalar)
        if strategy is None:
            sampler = self.tg.sample()
        else:
            sampler = self.tg.sample(strategy)
        sampler(rng, self.values, self.weights)
        assert numpy.isfinite(self.values).all()
        assert numpy.isfinite(self.weights).all()
        print(self.weights.min(), self.weights.max(), self.weights.sum())
        if False:
            self.ranges = numpy.array([self.values.min(axis=0),
                                       self.values.max(axis=0)]).transpose()
        else:
            self.ranges = numpy.array([[-1, 10], [-1, 10]], dtype=float)

    def eval1d(self, dim, x):
        i = self.dimensions.index(dim)
        y = numpy.linspace(0, self.ranges[not i, 1], 200)
        xg, yg = numpy.meshgrid(x, y)
        full = self.eval2d(dim, ('x' if dim == 'y' else 'y'), xg, yg)
        r = numpy.trapz(full, x=y, axis=0)
        return r

    def eval2d(self, xDim, yDim, x, y):
        evaluator = self.tg.evaluateLog()
        i = self.dimensions.index(xDim)
        j = self.dimensions.index(yDim)
        r = numpy.zeros(x.size, dtype=float)
        xy = numpy.zeros((x.size, 2), dtype=float)
        xy[:, i] = x.flat
        xy[:, j] = y.flat
        evaluator(xy, r)
        r -= self.tg.getLogIntegral()
        return numpy.exp(-r).reshape(x.shape)


def display(tg, nSamples=5000, strategy=None):
    data = TruncatedGaussianData(tg, nSamples=nSamples, strategy=strategy)
    figure = matplotlib.pyplot.figure(figsize=(10, 10))
    p = lsst.meas.modelfit.display.DensityPlot(figure, data)
    p.layers["hist"] = lsst.meas.modelfit.display.HistogramLayer()
    p.layers["samples"] = lsst.meas.modelfit.display.ScatterLayer()
    p.layers["func"] = lsst.meas.modelfit.display.SurfaceLayer()
    p.draw()
    return p

if __name__ == "__main__":
    if True:
        angle = -numpy.pi/6
        rot1 = lsst.afw.geom.LinearTransform.makeRotation(angle)
        rot2 = lsst.afw.geom.LinearTransform.makeRotation(-angle)
        scale = lsst.afw.geom.LinearTransform.makeScaling(1.0, 0)
        hessian = (rot2*scale*rot1).getMatrix()
        mu = numpy.array([2.0, 2.0])
        gradient = numpy.dot(hessian, -mu)
        tg = lsst.meas.modelfit.TruncatedGaussian.fromSeriesParameters(0.0, gradient, hessian)
    else:
        angle = numpy.pi/6
        rot1 = lsst.afw.geom.LinearTransform.makeRotation(angle)
        rot2 = lsst.afw.geom.LinearTransform.makeRotation(-angle)
        scale = lsst.afw.geom.LinearTransform.makeScaling(8.0, 2.0)
        sigma = (rot2*scale*rot1).getMatrix()
        mu = numpy.array([-3.0, 3.0])
        tg = lsst.meas.modelfit.TruncatedGaussian.fromStandardParameters(mu, sigma)
    p = display(tg)
    matplotlib.pyplot.show()
