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
"""
Config classes and registry for Bayesian priors
"""

import os
import numpy

import lsst.pex.config

from . import multifitLib

__all__ = ("priorRegistry", "registerPrior", "fitEllipticitySpline")

priorRegistry = lsst.pex.config.makeRegistry(
    """Registry for Bayesian priors on galaxy parameters

    The Configurables (callables that take a Config as their first argument)
    in the registry should return a subclass of the Prior class, and take one
    additional 'pixelScale' keyword argument (Angle/pixel).
    """
)

def registerPrior(name):
    """Decorator to add a Config class with a makePrior static method to the
    prior registry.
    """
    def decorate(cls):
        cls.makePrior.ConfigClass = cls
        priorRegistry.register(name, cls.makePrior)
        return cls
    return decorate

@registerPrior("flat")
class FlatPriorConfig(lsst.pex.config.Config):
    maxRadius = lsst.pex.config.Field(dtype=float, default=10.0, doc="Maximum radius in pixels")
    maxEllipticity = lsst.pex.config.Field(dtype=float, default=1.0,
                                           doc="Maximum ellipticity, in ReducedShear parametrization")

    @staticmethod
    def makePrior(config, pixelScale):
        return multifitLib.FlatPrior(config.maxRadius, config.maxEllipticity)

@registerPrior("mixture")
class MixturePriorConfig(lsst.pex.config.Config):
    filename = lsst.pex.config.Field(
        dtype=str, default="s13-v2-disk-08.fits",
        doc="Filename for mixture data file to load; relative to $MEAS_MULTIFIT_DIR/data unless absolute"
        )

    @staticmethod
    def makePrior(config, pixelScale):
        if os.path.isabs(config.filename):
            path = config.filename
        else:
            path = os.path.join(os.environ["MEAS_MULTIFIT_DIR"], "data", config.filename)
        mixture = multifitLib.Mixture3.readFits(path)
        # convert from log(r) in arcseconds to log(r) in pixels
        mixture.shift(2, -numpy.log(pixelScale.asArcseconds()))
        return multifitLib.MixturePrior(mixture)


def fitMixture(data, nComponents, minFactor=0.25, maxFactor=4.0, nIterations=20, df=float("inf")):
    """Fit a Mixture distribution to a set of (e1, e2, r) data points

    @param[in] data           array of data points to fit; shape=(N,3)
    @param[in] nComponents    number of components in the mixture distribution
    @param[in] minFactor      ellipticity variance of the smallest component in the initial mixture,
                              relative to the measured variance
    @param[in] maxFactor      ellipticity variance of the largest component in the initial mixture,
                              relative to the measured variance
    @param[in] nIterations    number of expectation-maximization update iterations
    @param[in] df             number of degrees of freedom for component Student's T distributions
                              (inf=Gaussian).
    """
    components = lsst.meas.multifit.Mixture3.ComponentList()
    rMu = data[:,2].mean()
    rSigma = data[:,2].var()
    eSigma = 0.5*(data[:,0].var() + data[:,1].var())
    mu = numpy.array([0.0, 0.0, rMu], dtype=float)
    baseSigma = numpy.array([[eSigma, 0.0, 0.0],
                             [0.0, eSigma, 0.0],
                             [0.0, 0.0, rSigma]])
    for factor in numpy.linspace(minFactor, maxFactor, nComponents):
        sigma = baseSigma.copy()
        sigma[:2,:2] *= factor
        components.append(lsst.meas.multifit.Mixture3.Component(1.0, mu, sigma))
    mixture = lsst.meas.multifit.Mixture3(components, df)
    restriction = lsst.meas.multifit.MixturePrior.getUpdateRestriction()
    for i in range(nIterations):
        mixture.updateEM(data, restriction)
    return mixture

def fitEllipticitySpline(ellipticity, histBins=100, outerKnotSpacing=0.01, nInteriorKnots=5, doPlot=False):
    import scipy.optimize
    import scipy.stats

    h, eg = numpy.histogram(ellipticity, bins=histBins, normed=True)
    de = (eg[1:] - eg[:-1]).mean()
    eh = 0.5*(eg[1:] + eg[:-1])
    knots = numpy.concatenate(
        ([0.0],
         numpy.linspace(outerKnotSpacing, eg[-1], nInteriorKnots),
         [eg[-1] + outerKnotSpacing, eg[-1] + 2*outerKnotSpacing]
         )
        )
    basis = multifitLib.ConstrainedSplineBasis(multifitLib.SplineBasis(knots, 3))
    basis.addConstraint(knots[0], 0.0, 1)  # constrain first derivative to zero at lower bound
    basis.addConstraint(knots[-1], 0.0, 1) # constrain first derivative to zero at upper bound
    basis.addConstraint(knots[-1], 0.0, 0) # constrain value to zero at upper bound
    basis.addIntegralConstraint()
    bh = numpy.zeros((eh.size, basis.getBasisSize()), dtype=float)
    ch = numpy.zeros((eh.size,), dtype=float)
    basis.evaluate(eh, bh, ch)

    h1 = h - h.min()
    kappa1 = (h1*de).sum() / (h*de).sum()
    def model1(e, x):
        return scipy.stats.t.pdf(e, x[1], scale=x[0])
    def func1(x):
        return kappa1*model1(eh, x) - h1
    sigma1, flags = scipy.optimize.leastsq(func1, numpy.array([0.3, 3.0]))

    eta1, residues, rank, sv = numpy.linalg.lstsq(bh, h.min() - ch)
    x1 = numpy.array([kappa1] + list(sigma1) + list(eta1))

    def modelh(x):
        return x[0] * model1(eh, x[1:3]) + (1.0 - x[0])*(numpy.dot(bh, x[3:]) + ch)
    def funch(x):
        return modelh(x) - h
    x2, flags = scipy.optimize.leastsq(funch, x1)

    sf = multifitLib.SplineFunction(basis.getSplineBasis(), basis.unconstrainCoefficients(x2[3:]))
    kappa = x2[0]
    tdist = scipy.stats.t(x2[2], scale=x2[1])

    if doPlot:
        from matplotlib import pyplot
        ep = numpy.linspace(0.0, 1.8, 5000)
        bp = numpy.zeros((ep.size, basis.getBasisSize()), dtype=float)
        cp = numpy.zeros((ep.size,), dtype=float)
        basis.evaluate(ep, bp, cp)
        def modelp(x):
            return x[0] * model1(ep, x[1:3]) + (1.0 - x[0])*(numpy.dot(bp, x[3:]) + cp)
        pyplot.plot(eh, h, 'ok')
        pyplot.plot(ep, kappa1*model1(ep, sigma1) + h.min(), '-r')
        pyplot.plot(ep, modelp(x1), '-b')
        pyplot.plot(ep, kappa*tdist.pdf(ep) + (1.0 - kappa)*sf(ep), '-m', linewidth=2)
        pyplot.plot(ep, modelp(x2), ':g')
        for k in knots:
            pyplot.axvline(k, color='c')
        pyplot.xlim(0, eg[-1]+0.2)
        pyplot.show()

    return kappa, tdisk, sf

