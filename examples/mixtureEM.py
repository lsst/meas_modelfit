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

import numpy
from matplotlib import pyplot
import mpl_toolkits.mplot3d

import lsst.meas.multifit

rng = lsst.afw.math.Random()

def makeRandomMixture(nDim, nComponents, df=float("inf")):
    l = lsst.meas.multifit.Mixture[nDim].ComponentList()
    for i in range(nComponents):
        mu = numpy.random.randn(nDim)*4
        a = numpy.random.randn(nDim+1,nDim)
        sigma = numpy.dot(a.transpose(),a) + numpy.identity(nDim)
        l.append(lsst.meas.multifit.Mixture[nDim].Component(numpy.random.rand(), mu, sigma))
    return lsst.meas.multifit.Mixture[nDim](l, df)

def initPlot1(fig, x, w):
    axes = fig.add_subplot(1,1,1)
    axes.hist(x, bins=100, normed=True, facecolor='k', alpha=0.2, range=(-15, 15), weights=w)
    xp = numpy.linspace(-15, 15, 2000)
    return axes, xp

def plotMixture1(initData, mixture, color, **kwds):
    axes, x = initData
    p = numpy.zeros(x.shape[0], dtype=float)
    mixture.evaluate(x.reshape(-1, 1), p)
    axes.plot(x, p, color, **kwds)

def finishPlot1(initData):
    axes, x = initData
    axes.legend()

def initPlot2(fig, x, w):
    axes2 = fig.add_subplot(2,2,1)
    h, xe, ye = numpy.histogram2d(x[:,0], x[:,1], weights=w, bins=(50, 50), normed=True,
                                  range=((-15, 15),(-15,15)))
    axes2.imshow(h.transpose(), origin='lower', interpolation='nearest',
                 extent=(-15, 15, -15, 15), cmap=pyplot.cm.Greys)
    xp = numpy.linspace(-15, 15, 2000)
    yp = numpy.linspace(-15, 15, 2000)
    gx, gy = numpy.meshgrid(xp, yp)
    gxy = numpy.array([gx.flatten(), gy.flatten()]).transpose().copy()
    axesY = fig.add_subplot(2,2,2, sharey=axes2)
    axesY.hist(x[:,1], bins=100, normed=True, range=(-15,15), facecolor='k', alpha=0.4, weights=w,
               orientation='horizontal', linewidth=0)
    axesX = fig.add_subplot(2,2,3, sharex=axes2)
    axesX.hist(x[:,0], bins=100, normed=True, range=(-15,15), facecolor='k', alpha=0.4, weights=w,
               linewidth=0)
    return axes2, axesY, axesX, xp, yp, gxy

def plotMixture2(initData, mixture, color, **kwds):
    axes2, axesY, axesX, xp, yp, gxy = initData
    pxy = numpy.zeros(gxy.shape[0], dtype=float)
    mixture.evaluate(gxy, pxy)
    axes2.contour(gxy[:,0].reshape(xp.size, yp.size), gxy[:,1].reshape(xp.size, yp.size),
                  -numpy.log10(pxy.reshape(xp.size, yp.size)), colors=color, alpha=0.6, **kwds)
    mixtureX = mixture.project(0)
    mixtureY = mixture.project(1)
    px = numpy.zeros(xp.size, dtype=float)
    py = numpy.zeros(yp.size, dtype=float)
    mixtureX.evaluate(xp.reshape(-1,1), px)
    mixtureY.evaluate(yp.reshape(-1,1), py)
    axesX.plot(xp, px, color, **kwds)
    axesY.plot(py, yp, color, **kwds)

def finishPlot2(initData):
    axes2, axesY, axesX, xp, yp, gxy = initData
    axesX.legend()
    axes2.set_ylim(-15, 15)
    axes2.set_xlim(-15, 15)

initPlot = {1: initPlot1, 2: initPlot2}
plotMixture = {1: plotMixture1, 2: plotMixture2}
finishPlot = {1: finishPlot1, 2: finishPlot2}

def doTestEM(inMixture, label, nSamples=100000, scatter=0.8, nIterations=50, importanceSample=False):
    nDim = inMixture.getDimension()
    outComponents = inMixture.ComponentList()
    for inComponent in inMixture:
        mu = (inComponent.getMu() + numpy.random.randn()*scatter).reshape(nDim)
        sigma = inComponent.getSigma() + numpy.identity(nDim)
        outComponents.append(inMixture.Component(1.0, mu, sigma))
    outMixture = type(inMixture)(outComponents, inMixture.getDegreesOfFreedom())

    x = numpy.zeros((nSamples, nDim), dtype=float)
    if importanceSample:
        outMixture.draw(rng, x)
        w = numpy.zeros(nSamples, dtype=float)
        inMixture.evaluate(x, w)
        q = numpy.zeros(nSamples, dtype=float)
        outMixture.evaluate(x, q)
        w /= q
    else:
        inMixture.draw(rng, x)
        w = numpy.ones(nSamples, dtype=float)
        w /= nSamples

    fig = pyplot.figure(label)
    plotData = initPlot[nDim](fig, x, w)
    plotMixture[nDim](plotData, inMixture, 'k', label='truth')
    plotMixture[nDim](plotData, outMixture, 'r', label='initial')
    for i in range(nIterations):
        outMixture.updateEM(x, w)
    plotMixture[nDim](plotData, outMixture, 'b', label='final')
    finishPlot[nDim](plotData)


if __name__ == "__main__":
    #m1 = makeRandomMixture(1, 3)
    #doTestEM(m1, "3x Gaussian, 1-d, uniform weights")
    #m2 = makeRandomMixture(1, 3, df=4)
    #doTestEM(m2, "3x Student's T, 1-d, uniform weights")
    m3 = makeRandomMixture(2, 3)
    doTestEM(m3, "3x Gaussian, 2-d, uniform weights")
    doTestEM(m3, "3x Gaussian, 2-d, importance weights", importanceSample=True)
    pyplot.show()
