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

__all__ = ("BaseSamplerConfig", "BaseSamplerTask", "ImportanceSamplerConfig", "AdaptiveImportanceSamplerConfig", "AdaptiveImportanceSamplerTask")

BaseSamplerConfig = lsst.pex.config.Config

class BaseSamplerTask(lsst.pipe.base.Task):
    """Base class for sampler subtask; responsible for constructing BaseSample subclass instances
    that do the actual work via their run() method.

    The same sampler can be used for both multi-epoch and single-epoch processing; the difference is
    abstracted away by the likelihood function object passed to run(), and to a lesser extent by the
    different initialization options provided by setup() and reset().
    """

    ConfigClass = BaseSamplerConfig

    def __init__(self, **kwds):
        lsst.pipe.base.Task.__init__(self, **kwds)

    def makeProposal(self, exposure, parameters):
        raise NotImplementedError("makeProposal not implemented for this sampler")

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
    initialSigma = lsst.pex.config.Field(
        dtype=float, default=0.4,
        doc="Initial width of proposal components"
        )
    initialSpacing = lsst.pex.config.Field(
        dtype=float, default=0.8,
        doc="Initial spacing of proposal components"
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

    def __init__(self, sampleSchema, **kwds):
        BaseSamplerTask.__init__(self, **kwds)
        self.rng = lsst.afw.math.Random(self.config.rngAlgorithm, self.config.rngSeed)
        self.impl = multifitLib.AdaptiveImportanceSampler(
            sampleSchema, self.rng, self.config.getIterationMap(), self.config.doSaveIterations
            )

    def makeLatinCube(self, nComponents, parameterDim):
        design = numpy.zeros((nComponents, parameterDim), dtype=float)
        numpy.random.seed(int(self.rng.uniformInt(1000)))
        x = numpy.linspace(-1, 1, nComponents)
        for j in xrange(parameterDim):
            design[:,j] = x[numpy.random.permutation(nComponents)]
        # note: we could do some permutations to make this a better sampling
        # of the space (i.e. reduce correlations, move things further apart),
        # but that gets complicated pretty fast, and so far it's simple and
        # easy
        return design

    def makeProposal(self, exposure, parameters):
        components = multifitLib.Mixture.ComponentList()
        sigma = numpy.identity(parameters.size, dtype=float) * self.config.initialSigma**2
        design = self.makeLatinCube(self.config.nComponents, parameters.size)
        for n in xrange(self.config.nComponents):
            mu = parameters.copy()
            mu[:] += design[n,:]*self.config.initialSpacing
            components.append(multifitLib.Mixture.Component(1.0, mu, sigma))
        df = self.config.degreesOfFreedom or float("inf")
        return multifitLib.Mixture(parameters.size, components, df)
