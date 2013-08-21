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

import lsst.pex.config
import lsst.pipe.base
import lsst.afw.geom.ellipses

from . import multifitLib

__all__ = ("BaseSamplerConfig", "BaseSamplerTask", "NaiveGridSamplerConfig", "NaiveGridSamplerTask")

BaseSamplerConfig = lsst.pex.config.Config

class BaseSamplerTask(lsst.pipe.base.Task):
    """Base class for sampler subtask; responsible for constructing BaseSample subclass instances
    that do the actual work via their run() method.

    The same sampler can be used for both multi-epoch and single-epoch processing; the difference is
    abstracted away by the objective function object passed to run(), and to a lesser extent by the
    different initialization options provided by setup() and reset().
    """

    ConfigClass = BaseSamplerConfig

    def __init__(self, **kwds):
        lsst.pipe.base.Task.__init__(self, **kwds)

    def setup(self, exposure, ellipse, center, prior):
        """Bootstrap the sampler for an object with no previous samples, using the given reference catalog
        ellipse (an afw.geom.ellipses.BaseCore object) and center position.

        @return an instance of a subclass of BaseSampler
        """
        raise NotImplementedError("setup() not implemented for this sampler")

    def reset(self, samples, center, prior):
        """Reinitialize the sampler state given a SampleSet from a previous sampler on the same object.

        The given SampleSet must have the same model definition as self.model.

        @return an instance of a subclass of BaseSampler
        """
        raise NotImplementedError("reset() not implemented for this sampler")

class NaiveGridSamplerConfig(BaseSamplerConfig):
    nRadiusSteps = lsst.pex.config.Field(
        dtype=int,
        default=12,
        doc="Number of radius steps in grid"
    )
    maxRadiusFactor = lsst.pex.config.Field(
        dtype=float,
        default=2.0,
        doc="Maximum radius in grid, as scaling factor multiplied by truth value radius"
    )
    maxEllipticity = lsst.pex.config.Field(
        dtype=float,
        default=0.9,
        doc="Maximum ellipticity magnitude in grid"
    )
    ellipticityStepSize = lsst.pex.config.Field(
        dtype=float,
        default=0.15,
        doc="ellipticity grid step size"
    )

class NaiveGridSamplerTask(BaseSamplerTask):
    ConfigClass = NaiveGridSamplerConfig

    def setup(self, exposure, ellipse, center, prior):
        axes = lsst.afw.geom.ellipses.Axes(ellipse)
        maxRadius = axes.getA() * self.config.maxRadiusFactor
        return multifitLib.NaiveGridSampler(
            center,
            self.config.nRadiusSteps,
            self.config.ellipticityStepSize,
            maxRadius,
            self.config.maxEllipticity
        )

@lsst.pex.config.wrap(multifitLib.ImportanceSamplerControl)
class ImportanceSamplerConfig(lsst.pex.config.Config):
    pass

class AdaptiveImportanceSamplerConfig(BaseSamplerConfig):
    rngAlgorithm = lsst.pex.config.Field(
        dtype=str, default="MT19937",
        doc="Algorithm used by pseudo-random number generator (see afw::math::Random)"
        )
    rngSeed = lsst.pex.config.Field(
        dtype=int, default=1,
        doc="Initial seed for pseudo-random number generator (see afw::math::Random)"
        )
    initialEllipticitySigma = lsst.pex.config.Field(
        dtype=float, default=0.4,
        doc="Initial width of proposal components in (ConformalShear) ellipticity dimensions"
        )
    initialEllipticitySpacing = lsst.pex.config.Field(
        dtype=float, default=0.8,
        doc="Initial spacing of proposal components in (ConformalShear) ellipticity dimensions"
        )
    initialRadiusSigma = lsst.pex.config.Field(
        dtype=float, default=0.5,
        doc="Initial width of proposal components in log(radius) dimensions"
        )
    initialRadiusSpacing = lsst.pex.config.Field(
        dtype=float, default=0.8,
        doc="Initial (random) scatter of proposal components in log(radius) dimensions"
        )
    nComponents = lsst.pex.config.Field(
        dtype=int, default=10, doc="Number of mixture components in proposal distribution"
        )
    degreesOfFreedom = lsst.pex.config.Field(
        dtype=float, optional=True, default=8.0,
        doc="Degrees-of-freedom for proposal Student's T distributions (None==inf==Gaussian)"
        )
    iterations = lsst.pex.config.ConfigDictField(
        keytype=int, itemtype=ImportanceSamplerConfig, default={},
        doc=("Sequence of importance sampling iterations, ordered by their keys, which should be")
        )
    doSaveIterations = lsst.pex.config.Field(
        dtype=bool, default=False,
        doc="Whether to save intermediate SampleSets and proposal distributions for debugging perposes"
        )

    def getIterationMap(self):
        """Transform the iterations config dict into a map of C++ control objects."""
        m = multifitLib.ImportanceSamplerControlMap()
        for k, v in self.iterations.items():
            m[k] = v.makeControl()
        return m

    def setDefaults(self):
        self.iterations[0] = ImportanceSamplerConfig(nUpdateSteps=1)
        self.iterations[1] = ImportanceSamplerConfig(nUpdateSteps=2, targetPerplexity=0.1, maxRepeat=2)
        self.iterations[2] = ImportanceSamplerConfig(nUpdateSteps=8, targetPerplexity=0.95, maxRepeat=3)

class AdaptiveImportanceSamplerTask(BaseSamplerTask):
    ConfigClass = AdaptiveImportanceSamplerConfig

    def __init__(self, **kwds):
        BaseSamplerTask.__init__(self, **kwds)
        self.rng = lsst.afw.math.Random(self.config.rngAlgorithm, self.config.rngSeed)

    def makeLatinCube(self, n):
        design = numpy.zeros((n,3), dtype=float)
        numpy.random.seed(int(self.rng.uniformInt(1000)))
        x = numpy.linspace(-1, 1, n)
        for j in xrange(3):
            design[:,j] = x[numpy.random.permutation(n)]
        # note: we could do some permutations to make this a better sampling
        # of the space (i.e. reduce correlations, move things further apart),
        # but that gets complicated pretty fast, and so far it's simple and
        # easy
        return design

    def setup(self, exposure, ellipse, center, prior):
        separable = lsst.afw.geom.ellipses.SeparableConformalShearLogTraceRadius(ellipse)
        fiducial = separable.getParameterVector()
        components = multifitLib.Mixture3.ComponentList()
        sigma = numpy.array([[self.config.initialEllipticitySigma**2, 0.0, 0.0],
                             [0.0, self.config.initialEllipticitySigma**2, 0.0],
                             [0.0, 0.0, self.config.initialRadiusSigma**2]], dtype=float)
        design = self.makeLatinCube(self.config.nComponents)
        for n in xrange(self.config.nComponents):
            mu = fiducial.copy()
            mu[:2] += design[n,:2]*self.config.initialEllipticitySpacing
            mu[2] += design[n,2]*self.config.initialRadiusSpacing
            components.append(multifitLib.Mixture3.Component(1.0, mu, sigma))
        df = self.config.degreesOfFreedom or float("inf")
        proposal = multifitLib.Mixture3(components, df)
        return multifitLib.AdaptiveImportanceSampler(self.rng, proposal, prior, center,
                                                     self.config.getIterationMap(),
                                                     self.config.doSaveIterations)

    def reset(self, samples, center, prior):
        proposal = samples.getProposal()
        return multifitLib.AdaptiveImportanceSampler(self.rng, proposal, prior, center,
                                                     self.config.getIterationMap(),
                                                     self.config.doSaveIterations)
