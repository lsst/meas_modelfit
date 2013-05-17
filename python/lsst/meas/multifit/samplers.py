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

import lsst.pex.config
import lsst.pipe.base
import lsst.afw.geom.ellipses

from .multifitLib import NaiveGridSampler

class BaseSamplerConfig(lsst.pipe.base.Task.ConfigClass):
    pass

class BaseSamplerTask(lsst.pipe.base.Task):
    """Base class for sampler subtask; responsible for constructing BaseSample subclass instances
    that do the actual work via their run() method.

    The same sampler can be used for both multi-epoch and single-epoch processing; the difference is
    abstracted away by the objective function object passed to run(), and to a lesser extent by the
    different initialization options provided by setup() and reset().
    """

    ConfigClass = BaseSamperConfig

    def __init__(self, schema, model, **kwds):
        lsst.pipe.base.Task.__init__(self, **kwds)
        self.model = model

    def setup(self, exposure, source):
        """Bootstrap the sampler for an object with no prior samples, using measurements from
        the given source and operations on the pixels of the given exposure.

        @return an instance of a subclass of BaseSampler
        """
        raise NotImplementedError("setup() not implemented for this sampler")

    def reset(self, samples):
        """Reinitialize the sampler state given a SampleSet from a previous sampler on the same object.

        The given SampleSet must have the same model definition as self.model.

        @return a lsst.pipe.base.Struct with:
          sampler: an instance of a subclass of BaseSampler
          footprint: an afw.detection.Footprint that defines the pixel region to fit
        """
        raise NotImplementedError("reset() not implemented for this sampler")

class NaiveGridSamplerConfig(BaseSamplerConfig):
    nRadiusSteps = lsst.pex.config.Field(
        dtype=int,
        default=20,
        doc="Number of radius steps in grid"
    )
    maxRadiusFactor = lsst.pex.config.Field(
        dtype=float,
        default=4.0,
        doc="Maximum radius in grid, as scaling factor multiplied by adaptive moments major axis ratio"
    )
    maxEllipticity = lsst.pex.config.Field(
        dtype=float,
        default=0.9,
        doc="Maximum ellipticity magnitude in grid"
    )
    ellipticityStepSize = lsst.pex.config.Field(
        dtype=float,
        default=0.5,
        doc="ellipticity grid step size"
    )

class NaiveGridSamplerTask(BaseSamplerTask):
    ConfigClass = NaiveGridSamplerConfig

    def setup(self, exposure, source):
        axes = lsst.afw.geom.ellipses.Axes(source.getShape())
        maxRadius = axes.getA() * self.config.maxRadiusFactor
        return NaiveGridSampler(
            source.getCentroid(),
            self.config.nRadiusSteps,
            self.config.ellipticityStepSize,
            maxRadius,
            self.config.maxEllipticity
        )
