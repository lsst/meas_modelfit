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

from . import modelfitLib

__all__ = ("ImportanceSamplerConfig", "AdaptiveImportanceSamplerConfig", "AdaptiveImportanceSamplerTask")

@lsst.pex.config.wrap(modelfitLib.ImportanceSamplerControl)
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
    maxRetries = lsst.pex.config.Field(
        dtype=int, default=0,
        doc="Number of times we attempt to fit an object with different RNG states before giving up"
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
        m = modelfitLib.ImportanceSamplerControlMap()
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

    def __init__(self, schema, keys, model, prior, previous=None, **kwds):
        lsst.pipe.base.Task.__init__(self, **kwds)
        # n.b. schema argument is for modelfits catalog; self.sampleSchema is for sample catalog
        self.sampleSchema = lsst.afw.table.Schema()
        self.rng = lsst.afw.math.Random(self.config.rngAlgorithm, self.config.rngSeed)
        if self.config.doMarginalizeAmplitudes:
            self.interpreter = modelfitLib.MarginalSamplingInterpreter(self.sampleSchema, model, prior)
        else:
            self.interpreter = modelfitLib.DirectSamplingInterpreter(self.sampleSchema, model, prior)
        self.sampler = modelfitLib.AdaptiveImportanceSampler(
            self.sampleSchema, self.rng, self.config.getIterationMap(), self.config.doSaveIterations
            )
        self.keys = keys
        self.keys["rngstate"] = schema.addField(
            "rngstate", type=str, size=self.rng.getStateSize(),
            doc="Random number generator state (blob)"
            )
        self.previous = previous
        if self.previous is None: return
        if not isinstance(self.previous, AdaptiveImportanceSamplerTask):
            raise NotImplementedError("Cannot use previous catalog with"
                                      " non-AdaptiveImportanceSampler fitter")
        if not self.previous.config.doMarginalizeAmplitudes and self.config.doMarginalizeAmplitudes:
            raise NotImplementedError("Cannot use previous catalog with direct sampling to warm-start"
                                      " marginal sampling")
        if self.previous.config.doMarginalizeAmplitudes and not self.config.doMarginalizeAmplitudes:
            self.unnest = modelfitLib.UnnestMarginalSamples(
                self.previous.sampleSchema,
                self.sampleSchema,
                self.previous.interpreter,
                self.interpreter,
                self.rng
                )
        else:
            self.unnest = None

    def makeSampleTable(self):
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

    def initialize(self, outRecord):
        """Initialize an output record, setting any derived fields and record
        attributes (i.e. samples or pdf) needed before calling run().

        This method is not called when using a "warm start" from a previous fit.
        """
        parameters = numpy.zeros(self.interpreter.getParameterDim(), dtype=modelfitLib.Scalar)
        self.interpreter.packParameters(outRecord[self.keys['initial.nonlinear']],
                                        outRecord[self.keys['initial.amplitudes']],
                                        parameters)
        components = modelfitLib.Mixture.ComponentList()
        sigma = numpy.identity(parameters.size, dtype=float) * self.config.initialSigma**2
        design = self.makeLatinCube(self.rng, self.config.nComponents, parameters.size)
        for n in xrange(self.config.nComponents):
            mu = parameters.copy()
            mu[:] += design[n,:]*self.config.initialSpacing
            components.append(modelfitLib.Mixture.Component(1.0, mu, sigma))
        df = self.config.degreesOfFreedom or float("inf")
        proposal = modelfitLib.Mixture(parameters.size, components, df)
        outRecord.setPdf(proposal)

    def adaptPrevious(self, prevRecord, outRecord):
        """Adapt a previous record (fit using self.previous as the fitter task), filling in the
        fields and attributes of outRecord to put it in a state ready for run().
        """
        if self.unnest is not None:
            self.unnest.apply(prevRecord, outRecord)

    def run(self, likelihood, outRecord):
        """Do the actual fitting, using the given likelihood, update the 'pdf' and 'samples' attributes,
        and save best-fit values in the 'fit.parameters' field.
        """
        objective = modelfitLib.makeSamplingObjective(self.interpreter, likelihood)
        nRetries = 0
        while True:
            try:
                outRecord.setString(self.keys["rngstate"], self.rng.getState())
                self.sampler.run(objective, outRecord.getPdf(), outRecord.getSamples())
                break
            except Exception as err:
                if nRetries >= self.config.maxRetries:
                    raise
                self.log.warn("Failure fitting object %s; retrying with new RNG state" % outRecord.getId())
                self.initialize(outRecord)
                nRetries += 1
