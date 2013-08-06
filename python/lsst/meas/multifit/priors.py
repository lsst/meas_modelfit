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

__all__ = ("priorRegistry", "registerPrior")

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

