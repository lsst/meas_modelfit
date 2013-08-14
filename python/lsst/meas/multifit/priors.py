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

import numpy

import lsst.pex.config

from . import multifitLib

__all__ = ("priorRegistry", "registerPrior")

priorRegistry = lsst.pex.config.makeRegistry(
    """Registry for Bayesian priors on galaxy parameters

    The Configurables (callables that take a Config as their first and, in this case, only argument)
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

@registerPrior("flat")
class FlatPriorConfig(lsst.pex.config.Config):
    maxRadius = lsst.pex.config.Field(dtype=float, default=10.0, doc="Maximum radius in pixels")
    maxEllipticity = lsst.pex.config.Field(dtype=float, default=1.0,
                                           doc="Maximum ellipticity, in ReducedShear parametrization")

    @staticmethod
    def makePrior(config):
        return multifitLib.FlatPrior(config.maxRadius, config.maxEllipticity)
