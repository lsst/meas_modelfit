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

from . import modelfitLib

__all__ = ("priorRegistry", "registerPrior")

priorRegistry = lsst.pex.config.makeRegistry(
    """Registry for Bayesian priors on galaxy parameters

    The Configurables (callables that take a Config as their first argument)
    in the registry should return a subclass of the Prior class.
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


@registerPrior("mixture")
class MixturePriorConfig(lsst.pex.config.Config):
    filename = lsst.pex.config.Field(
        dtype=str, default="s13-v2-disk-08.fits",
        doc="Filename for mixture data file to load; relative to $MEAS_MODELFIT_DIR/data unless absolute"
    )

    @staticmethod
    def makePrior(config):
        if os.path.isabs(config.filename):
            path = config.filename
        else:
            path = os.path.join(os.environ["MEAS_MODELFIT_DIR"], "data", config.filename)
        mixture = modelfitLib.Mixture.readFits(path)
        return modelfitLib.MixturePrior(mixture, "single-ellipse")


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
    components = lsst.meas.modelfit.Mixture.ComponentList()
    rMu = data[:, 2].mean()
    rSigma = data[:, 2].var()
    eSigma = 0.5*(data[:, 0].var() + data[:, 1].var())
    mu = numpy.array([0.0, 0.0, rMu], dtype=float)
    baseSigma = numpy.array([[eSigma, 0.0, 0.0],
                             [0.0, eSigma, 0.0],
                             [0.0, 0.0, rSigma]])
    for factor in numpy.linspace(minFactor, maxFactor, nComponents):
        sigma = baseSigma.copy()
        sigma[:2, :2] *= factor
        components.append(lsst.meas.modelfit.Mixture.Component(1.0, mu, sigma))
    mixture = lsst.meas.modelfit.Mixture(3, components, df)
    restriction = lsst.meas.modelfit.MixturePrior.getUpdateRestriction()
    for i in range(nIterations):
        mixture.updateEM(data, restriction)
    return mixture
