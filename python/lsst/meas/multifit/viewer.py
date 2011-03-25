#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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

import lsst.afw.image
import lsst.afw.detection
import lsst.afw.image
import lsst.afw.geom.ellipses
import lsst.meas.multifit
import numpy

from matplotlib import pyplot

class ViewerBase(object):

    def __init__(self, evaluator, footprint, parameters=None, coefficients=None):
        self.evaluator = evaluator
        self.footprint = footprint
        self.parameters = numpy.zeros(self.evaluator.getParameterSize(), dtype=float)
        self.coefficients = numpy.zeros(self.evaluator.getCoefficientSize(), dtype=float)
        self.modelMatrix = numpy.zeros((self.footprint.getArea(), self.evaluator.getCoefficientSize()),
                                       dtype=float)
        self.bbox = self.footprint.getBBox()
        self.dataImage = lsst.afw.image.ImageD(self.bbox)
        self.modelImage = lsst.afw.image.ImageD(self.bbox)
        self.residualImage = lsst.afw.image.ImageD(self.bbox)
        self.dataVector = self.evaluator.getDataVector()
        self.vmin = self.dataVector.min()
        self.vmax = self.dataVector.max()
        lsst.afw.detection.expandArray(
            self.footprint, self.dataVector, self.dataImage.getArray(), self.bbox.getMin()
            )
        if parameters is None:
            parameters = numpy.zeros(self.evaluator.getParameterSize(), dtype=float)
            self.evaluator.writeInitialParameters(parameters)
        self.update(parameters, coefficients)

    def update(self, parameters, coefficients=None):
        self.parameters[:] = parameters
        self.evaluator.evaluateModelMatrix(self.modelMatrix, self.parameters)
        r = None
        if coefficients is not None:
            self.coefficients[:] = coefficients
        else:
            r = self.solve()
        self.modelVector = numpy.dot(self.modelMatrix, self.coefficients)
        self.residualVector = self.dataVector - self.modelVector
        lsst.afw.detection.expandArray(
            self.footprint, self.modelVector, self.modelImage.getArray(), self.bbox.getMin()
            )
        lsst.afw.detection.expandArray(
            self.footprint, self.residualVector, self.residualImage.getArray(), self.bbox.getMin()
            )
        return r

    def solve(self):
        x, residues, rank, sv = numpy.linalg.lstsq(self.modelMatrix, self.evaluator.getDataVector())
        self.coefficients[:] = x
        return 0.5 * residues[0] + 0.5 * numpy.log(sv.sum())

    def plot(self, fignum=None):
        dbox = lsst.afw.geom.Box2D(self.bbox)
        extent = (dbox.getMinX(), dbox.getMaxX(), dbox.getMinY(), dbox.getMaxY())
        if fignum is None:
            figure = pyplot.figure()
        else:
            figure = pyplot.figure(fignum)
        def doCell(image, title):            
            pyplot.imshow(image.getArray(), origin='lower', interpolation='nearest',
                          extent=extent, vmin=self.vmin, vmax=self.vmax)
            self.plotGeometry()
            pyplot.title(title)
            pyplot.xlim(extent[0], extent[1])
            pyplot.ylim(extent[2], extent[3])
        pyplot.subplot(1, 3, 1)
        doCell(self.dataImage, "DATA")
        modelAxes = pyplot.subplot(1, 3, 2)
        doCell(self.modelImage, "MODEL")
        residualAxes = pyplot.subplot(1, 3, 3)
        doCell(self.residualImage, "RESIDUAL")
        cax = figure.add_axes([0.05, 0.08, 0.90, 0.04])
        pyplot.colorbar(cax=cax, orientation='horizontal')
        
class StarViewer(ViewerBase):

    @staticmethod
    def makeExample(sn=10.0):
        bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(0, 0), lsst.afw.geom.Extent2I(50, 60))
        footprint = lsst.afw.detection.Footprint(bbox)
        point = lsst.afw.geom.Point2D(25.0, 30.0)
        psf = lsst.afw.detection.createPsf("DoubleGaussian", 19, 19, 2.0, 1.0)
        localPsf = psf.getLocalPsf(point)
        vector = numpy.zeros(footprint.getArea(), dtype=float)
        localPsf.evaluatePointSource(footprint, vector, lsst.afw.geom.Extent2D())
        signal = vector.max()
        sigma = signal / sn
        exposure = lsst.afw.image.ExposureD(bbox)
        lsst.afw.detection.expandArray(
            footprint, vector, exposure.getMaskedImage().getImage().getArray(), bbox.getMin()
            )
        exposure.getMaskedImage().getImage().getArray()[:,:] \
            += numpy.random.normal(scale=sigma, size=(bbox.getHeight(), bbox.getWidth()))
        exposure.getMaskedImage().getVariance().getArray()[:,:] = sigma
        exposure.setPsf(psf)
        evaluator = lsst.meas.multifit.Evaluator.make(exposure, footprint, point, False, False)
        return StarViewer(evaluator, footprint)

    def update(self, parameters, coefficients=None):
        r = ViewerBase.update(self, parameters, coefficients)
        self.point = self.evaluator.extractPoint(0, self.parameters)
        return r

    def plotGeometry(self):
        pyplot.plot([self.point.getX()], [self.point.getY()], 'kx')
        
class GalaxyViewer(ViewerBase):

    def update(self, parameters, coefficients=None):
        r = ViewerBase.update(self, parameters, coefficients)
        self.ellipse = self.evaluator.extractEllipse(0, self.parameters)
        return r

    def plotGeometry(self):
        self.ellipse.plot(fill=False)

def makeGaussianEllipse(component, i, j):
    mu = component.getMu()
    sigma = component.getSigma()
    covariance = numpy.dot(sigma, sigma.transpose())
    print covariance
    return lsst.afw.geom.ellipses.Ellipse(
        lsst.afw.geom.ellipses.Quadrupole(float(covariance[j,j]), float(covariance[i,i]), 
                                          float(covariance[i,j])),
        lsst.afw.geom.Point2D(float(mu[j]), float(mu[i]))
        )

def plotSamples(table, importance):
    nParameters = table["parameters"].shape[1]
    mid = lambda x: 0.5*(x[:-1] + x[1:])
    for i in range(nParameters):
        for j in range(i + 1):
            pyplot.subplot(nParameters, nParameters, i * nParameters + j + 1)
            if i == j:
                h, e = numpy.histogram(table["parameters"][:,i], bins=20, weights=table["weight"], new=True)
                pyplot.plot(mid(e), h)
            else:
                h, xe, ye = numpy.histogram2d(table["parameters"][:,j], table["parameters"][:,i],
                                              bins=10, weights=table["weight"])
                pyplot.pcolor(xe, ye, h)
                #pyplot.plot(table["parameters"][:,j], table["parameters"][:,i], 'k,')
                for component in importance.getComponents():
                    ellipse = makeGaussianEllipse(component, i, j)
                    ellipse.plot(fill=False)
