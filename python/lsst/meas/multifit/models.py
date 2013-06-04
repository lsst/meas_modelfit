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
Config classes used to define various galaxy models.
"""

import numpy

import lsst.pex.config
import lsst.shapelet.tractor

__all__ = ("modelRegistry", "registerModel")

modelRegistry = lsst.pex.config.makeRegistry(
    """Registry for galaxy model definitions

    A galaxy definition is a Configurable (a callable that takes a Config object
    as its first and, in this case, only argument) that returns a
    shapelet.MultiShapeletBasis.
    """
)

def registerModel(name):
    """Decorator to add a Config class with a makeBasis static method to the
    model registry.
    """
    def decorate(cls):
        cls.makeBasis.ConfigClass = cls
        registry.register(name, cls.makeBasis)
        return cls
    return decorate

@registerModel("gaussian")
class GaussianModelConfig(lsst.pex.config.Config):
    """Config class used to define a simple Gaussian galaxy model.
    """
    radius = lsst.pex.config.Field(
        "Radius parameter to use for the model, in units of the RMS size (sigma)",
        dtype=float, default=1.0,
        )

    @staticmethod
    def makeBasis(config):
        """Create and return a MultiShapeletBasis corresponding to the config."""
        basis = lsst.shapelet.MultiShapeletBasis(1)
        basis.addComponent(config.radius, 0, numpy.array([[1.0]], dtypef=float))
        return basis

@registerModel("tractor")
class TractorModelConfig(lsst.pex.config.Config):
    """Config class used to define a MultiShapeletBasis approximation to a Sersic or Sersic-like profile,
    as optimized by Hogg and Lang's The Tractor.

    See the lsst.shapelet.tractor module for more information.
    """
    profile = lsst.pex.config.Field(
        "name of the profile: one of 'exp', 'dev', 'ser2', 'ser3', 'luv', or 'lux'",
        dtype=str, default='lux',
        )
    nComponents = lsst.pex.config.Field(
        "number of Gaussian components in the approximation",
        dtype=int, default=8,
        )
    maxRadius = lsst.pex.config.Field(
        ("maximum radius for which the multi-Gaussian approximation was optimized, in units"
         " of the half-light radius; None will use the profile-dependent default"),
        dtype=int, optional=True,
        )

    @staticmethod
    def makeBasis(config):
        """Load and return a MultiShapeletBasis corresponding to the config."""
        return lsst.shapelet.tractor.loadBasis(
            profile=config.profile,
            nComponents=config.nComponents,
            maxRadius=config.maxRadius
            )

@registerModel("bulge+disk")
class BulgeDiskModelConfig(lsst.pex.config.Config):
    """Config that defines the model used in two-component galaxy model fits.
    """
    disk = lsst.pex.config.ConfigurableField(
        "multi-Gaussian approximation to be used for the disk component of the model",
        target=TractorModelConfig.makeBasis,
        )
    bulge = lsst.pex.config.ConfigurableField(
        "multi-Gaussian approximation to be used for the bulge component of the model",
        target=TractorModelConfig.makeBasis,
        )
    bulgeRadius = lsst.pex.config.Field(
        "half-light radius of bulge in units of half-light radius of disk",
        dtype=float,
        default=0.6,
        )

    def setDefaults(self):
        self.disk.profile = "lux"
        self.bulge.profile = "luv"

    @staticmethod
    def makeBasis(config):
        """Return a MultiShapeletBasis with both disk and bulge components.
        """
        bulge = config.bulge.apply()
        bulge.scale(config.bulgeRadius)
        disk = config.disk.apply()
        disk.merge(bulge)
        return disk
