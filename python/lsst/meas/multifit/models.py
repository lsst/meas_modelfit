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

from . import multifitLib

__all__ = ("modelRegistry", "registerModel")

modelRegistry = lsst.pex.config.makeRegistry(
    """Registry for galaxy model definitions

    A galaxy model definition is a Configurable (a callable that takes a Config object
    as its first argument) that returns a lsst.meas.multifit.Model.  It also takes
    the center position of the source as a second argument.
    """
)

def registerModel(name):
    """Decorator to add a Config class with a makeModel static method to the
    model registry.
    """
    def decorate(cls):
        cls.makeModel.ConfigClass = cls
        modelRegistry.register(name, cls.makeModel)
        return cls
    return decorate

def getCenterEnum(config):
    """Helper function to turn config boolean option into enum value"""
    return multifitLib.Model.FIXED_CENTER if config.fixCenter else multifitLib.Model.SINGLE_CENTER

@registerModel("gaussian")
class GaussianModelConfig(lsst.pex.config.Config):
    """Config class used to define a simple Gaussian galaxy model.
    """
    radius = lsst.pex.config.Field(
        "Radius parameter to use for the model, in units of the RMS size (sigma)",
        dtype=float, default=1.0,
        )
    fixCenter = lsst.pex.config.Field(
        "Fix the center to the position derived from a previous centeroider?",
        dtype=bool, default=True
        )

    @staticmethod
    def makeModel(config):
        return multifitLib.Model.makeGaussian(getCenterEnum(config), self.config.radius)

class FixedSersicConfig(lsst.pex.config.Config):
    """Config class used to define a MultiShapeletBasis approximation to a Sersic or Sersic-like profile,
    as optimized by Hogg and Lang's The Tractor.  Intended for use as a subclass or nested config only,
    not a top-level model config.

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

    def makeBasis(self):
        """Return a MultiShapeletBasis corresponding to the config."""
        maxRadius = self.maxRadius if self.maxRadius is not None else 0
        return lsst.shapelet.RadialProfile.get(self.profile).getBasis(self.nComponents, maxRadius)

@registerModel("fixed-sersic")
class FixedSersicModelConfig(FixedSersicConfig):
    """A single-component fixed-index Sersic model, using a multi-Gaussian approximation to
    the profile (see FixedSersicConfig).
    """
    fixCenter = lsst.pex.config.Field(
        "Fix the center to the position derived from a previous centeroider?",
        dtype=bool, default=True
        )

    @staticmethod
    def makeModel(config):
        return multifitLib.Model.make(config.makeBasis(), getCenterEnum(config))

@registerModel("bulge+disk")
class BulgeDiskModelConfig(lsst.pex.config.Config):
    """Config that defines the model used in two-component galaxy model fits.
    """
    disk = lsst.pex.config.ConfigField(
        "multi-Gaussian approximation to be used for the disk component of the model",
        dtype=FixedSersicConfig
        )
    bulge = lsst.pex.config.ConfigField(
        "multi-Gaussian approximation to be used for the bulge component of the model",
        dtype=FixedSersicConfig
        )
    bulgeRadius = lsst.pex.config.Field(
        ("Half-light radius of bulge in units of half-light radius of disk. "
         "If None, the two components will have completely independent radii "
         "and ellipticities."),
        dtype=float,
        default=0.6,
        optional=True
        )
    fixCenter = lsst.pex.config.Field(
        "Fix the center to the position derived from a previous centeroider?",
        dtype=bool, default=True
        )

    def setDefaults(self):
        self.disk.profile = "lux"
        self.bulge.profile = "luv"

    @staticmethod
    def makeModel(config):
        bulge = config.bulge.makeBasis()
        disk = config.disk.makeBasis()
        if config.bulgeRadius is None:
            basisVector = multifitLib.Model.BasisVector()
            basisVector.append(disk)
            basisVector.append(bulge)
            prefixes = multifitLib.Model.NameVector()
            prefixes.append("exp.")
            prefixes.append("dev.")
            return multifitLib.Model.make(basisVector, prefixes, getCenterEnum(config))
        else:
            basis = disk
            bulge.scale(config.bulgeRadius)
            basis.merge(bulge)
            return multifitLib.Model.make(basis, getCenterEnum(config))
