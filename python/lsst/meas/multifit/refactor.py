
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

    ConfigClass = AdaptiveImportanceSamplerConfig

    def __init__(self, keys, model, prior, **kwds):
        # n.b. keys argument is for modelfits catalog; self.schema is for sample catalog
        self.schema = lsst.afw.table.Schema()
        self.rng = lsst.afw.math.Random(self.config.rngAlgorithm, self.config.rngSeed)
        self.objectiveFactory = multifitLib.SamplerObjectiveFactory(
            self.schema, model, prior, self.config.doMarginalizeAmplitudes
            )
        self.sampler = multifitLib.AdaptiveImportanceSampler(
            self.schema, self.rng, self.config.getIterationMap(), self.config.doSaveIterations
            )
        self.keys = keys

    @staticmethod
    def makeLatinCube(rng, nComponents, parameterDim):
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

    def initialize(self, model, record):
        parameters = record[self.keys["ref.parameters"]]
        self.objectiveFactory.mapParameters(record[self.keys["ref.nonlinear"]],
                                            record[self.keys["ref.amplitudes"]],
                                            parameters)
        components = multifitLib.Mixture.ComponentList()
        sigma = numpy.identity(parameters.size, dtype=float) * self.config.initialSigma**2
        design = self.makeLatinCube(self.rng, self.config.nComponents, model.getParameterDim())
        for n in xrange(self.config.nComponents):
            mu = parameters.copy()
            mu[:] += design[n,:]*self.config.initialSpacing
            components.append(multifitLib.Mixture.Component(1.0, mu, sigma))
        df = self.config.degreesOfFreedom or float("inf")
        proposal = multifitLib.Mixture(model.getParameterDim(), components, df)
        record.setPdf(proposal)

    def run(self, likelihood, record):
        objective = self.objectiveFactory(likelihood)
        self.sampler.run(objective, record.getPdf(), record.getSamples())
        # TODO: compute and set best-fit parameters
