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
import lsst.afw.table
import lsst.afw.math

from . import multifitLib

__all__ = ("ImportanceSamplerConfig", "AdaptiveImportanceSamplerConfig", "AdaptiveImportanceSamplerTask")

@lsst.pex.config.wrap(multifitLib.ImportanceSamplerControl)
class ImportanceSamplerConfig(lsst.pex.config.Config):
    pass

class AdaptiveImportanceSamplerConfig(lsst.pex.config.Config):
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
    doMarginalizeAmplitudes = lsst.pex.config.Field(
        dtype=bool,
        default=True,
        doc="Marginalize over amplitudes numerically instead of sampling them?"
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

class AdaptiveImportanceSamplerTask(lsst.pipe.base.Task):
    """A 'fitter' subtask for Measure tasks that uses adaptive importance sampling.
    """

    ConfigClass = AdaptiveImportanceSamplerConfig

    def __init__(self, schema, keys, model, prior, **kwds):
        lsst.pipe.base.Task.__init__(self, **kwds)
        # n.b. schema argument is for modelfits catalog; self.sampleSchema is for sample catalog
        self.sampleSchema = lsst.afw.table.Schema()
        self.rng = lsst.afw.math.Random(self.config.rngAlgorithm, self.config.rngSeed)
        if self.config.doMarginalizeAmplitudes:
            self.interpreter = multifitLib.MarginalSamplingInterpreter(self.sampleSchema, model, prior)
        else:
            self.interpreter = multifitLib.DirectSamplingInterpreter(self.sampleSchema, model, prior)
        self.sampler = multifitLib.AdaptiveImportanceSampler(
            self.sampleSchema, self.rng, self.config.getIterationMap(), self.config.doSaveIterations
            )
        self.keys = keys
        self.keys["ref.parameters"] = schema.addField(
            "ref.parameters", type="ArrayD", size=self.interpreter.getParameterDim(),
            doc="sampler parameters from reference catalog"
            )
        self.keys["fit.parameters"] = schema.addField(
            "fit.parameters", type="ArrayD", size=self.interpreter.getParameterDim(),
            doc="best-fit sampler parameters"
            )
        self.keys["rngstate"] = schema.addField(
            "rngstate", type=str, size=self.rng.getStateSize(),
            doc="Random number generator state (blob)"
            )

    def makeTable(self):
        """Return a Table object that can be used to construct sample records.
        """
        return lsst.afw.table.BaseTable.make(self.sampleSchema)

    @staticmethod
    def makeLatinCube(rng, nComponents, parameterDim):
        """Create a (nComponents, parameterDim) design matrix that represents
        an Latin hypercube sampling of a space with dimension parameterDim
        using nComponents samples.  The range of each dimension is assumed to
        be (-1,1).
        """
        design = numpy.zeros((nComponents, parameterDim), dtype=float)
        numpy.random.seed(int(rng.uniformInt(1000)))
        x = numpy.linspace(-1, 1, nComponents)
        for j in xrange(parameterDim):
            design[:,j] = x[numpy.random.permutation(nComponents)]
        # note: we could do some permutations to make this a better sampling
        # of the space (i.e. reduce correlations, move things further apart),
        # but that gets complicated pretty fast, and so far it's simple and
        # easy
        return design

    def initialize(self, record):
        """Initialize an output record, setting any derived fields and record
        attributes (i.e. samples or pdf) needed before calling run().

        This method is not called when using a "warm start" from a previous fit.
        """
        parameters = record[self.keys["ref.parameters"]]
        self.interpreter.packParameters(record[self.keys["ref.nonlinear"]],
                                        record[self.keys["ref.amplitudes"]],
                                        parameters)
        components = multifitLib.Mixture.ComponentList()
        sigma = numpy.identity(parameters.size, dtype=float) * self.config.initialSigma**2
        design = self.makeLatinCube(self.rng, self.config.nComponents, parameters.size)
        for n in xrange(self.config.nComponents):
            mu = parameters.copy()
            mu[:] += design[n,:]*self.config.initialSpacing
            components.append(multifitLib.Mixture.Component(1.0, mu, sigma))
        df = self.config.degreesOfFreedom or float("inf")
        proposal = multifitLib.Mixture(parameters.size, components, df)
        record.setPdf(proposal)

    def run(self, likelihood, record):
        """Do the actual fitting, using the given likelihood, update the 'pdf' and 'samples' attributes,
        and save best-fit values in the 'fit.parameters' field.
        """
        objective = multifitLib.makeSamplingObjective(self.interpreter, likelihood)
        record.setString(self.keys["rngstate"], self.rng.getState())
        self.sampler.run(objective, record.getPdf(), record.getSamples())
        record[self.keys["fit.parameters"]] = self.interpreter.computeParameterMean(record)
