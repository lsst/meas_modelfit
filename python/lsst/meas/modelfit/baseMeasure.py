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
import lsst.afw.image
import lsst.afw.geom
import lsst.afw.table

from . import modelfitLib
from .models import modelRegistry
from .priors import priorRegistry
from .optimizer import OptimizerTask
from .fitRegion import setupFitRegion

__all__ = ("BaseMeasureConfig", "BaseMeasureTask")

class BaseMeasureConfig(lsst.pex.config.Config):
    fitter = lsst.pex.config.ConfigurableField(
        target=OptimizerTask,
        doc="Subtask that actually does the fitting"
    )
    model = modelRegistry.makeField(
        default="fixed-sersic",
        doc="Definition of the galaxy model to fit"
    )
    prior = priorRegistry.makeField(
        default="mixture",
        doc="Bayesian prior on galaxy parameters"
    )
    fitRegion = lsst.pex.config.ConfigField(
        dtype=setupFitRegion.ConfigClass,
        doc="Parameters that control which pixels to include in the model fit"
    )
    progressChunk = lsst.pex.config.Field(
        dtype=int,
        default=100,
        doc="Show progress log message every [progressChunk] objects"
    )
    psf = lsst.pex.config.ConfigField(
        dtype=modelfitLib.PsfFitterConfig,
        doc="Config options for approximating the PSF using shapelets"
    )
    prepOnly = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        doc="If True, only prepare the catalog (match, transfer fields, fit PSF)"
    )
    maxObjects = lsst.pex.config.Field(
        dtype=int,
        default=None,
        optional=True,
        doc="If not None, clip the catalog and process only this many objects (for fast-debug purposes)"
    )
    doRaise = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        doc="If True, raise exceptions when individual objects fail instead of warning."
    )
    tag = lsst.pex.config.Field(
        dtype=str,
        default="none",
        doc=("A user-defined string added to the data ID of outputs, used to version different"
             " configurations in a way that allows them to be chained together.")
    )
    previous = lsst.pex.config.Field(
        dtype=str,
        default=None,
        optional=True,
        doc="The tag of a previous run of a meas_modelfit measure task to use as inputs"
    )

class BaseMeasureTask(lsst.pipe.base.CmdLineTask):
    """An intermediate base class for top-level model-fitting tasks.

    Subclasses (for different kinds of input data) must implement three methods:
    readInputs(), prepCatalog(), makeLikelihood(), and writeOutputs().

    Different fitting/sampling algorithms can be plugged in via the fitter subtask.
    """

    RunnerClass = lsst.pipe.base.ButlerInitializedTaskRunner
    ConfigClass = BaseMeasureConfig

    def __init__(self, butler=None, **kwds):
        """Initialize the measurement task, including the modelfits catalog schema,
        the model, prior, and calib objects, and the fitter subtask.
        """
        lsst.pipe.base.CmdLineTask.__init__(self, **kwds)
        # Because we're doing all the fitting in celestial coordinates (actually local tangent planes
        # that's aligned with celestial coordinates), we can define the model and prior
        # up front, and use them without modification for all Objects
        self.model = self.config.model.apply()
        self.prior = self.config.prior.apply()
        # now we set up the schema; this will be the same regardless of whether or not
        # we do a warm start (use a previous modelfits catalog for initial values).
        self.schema = modelfitLib.ModelFitTable.makeMinimalSchema()
        self.keys = {}
        self.keys["center"] = self.schema.addField(
            "center", type="PointD",
            doc="input centroid of the object in image coordinates"
            )
        self.keys["initial.nonlinear"] = self.schema.addField(
            "initial.nonlinear", type="ArrayD", size=self.model.getNonlinearDim(),
            doc="initial (pre-fit) nonlinear parameters"
            )
        self.keys["initial.amplitudes"] = self.schema.addField(
            "initial.amplitudes", type="ArrayD", size=self.model.getAmplitudeDim(),
            doc="initial (pre-fit) linear amplitudes"
            )
        self.keys["fixed"] = self.schema.addField(
            "fixed", type="ArrayD", size=self.model.getFixedDim(),
            doc="fixed nonlinear parameters"
            )
        self.keys["snr"] = self.schema.addField(
            "snr", type=float,
            doc="signal to noise ratio from source apFlux/apFluxErr"
            )
        self.keys["sys.position"] = self.schema.addField(
            "sys.position", type="Coord",
            doc="nominal position used to construct UnitSystem for parameters"
            )
        self.keys["sys.magnitude"] = self.schema.addField(
            "sys.mag", type=float,
            doc="nominal magnitude used to construct UnitSystem for parameters"
            )
        self.keys["fit.nonlinear"] = self.schema.addField(
            "fit.nonlinear", type="ArrayD", size=self.model.getNonlinearDim(),
            doc="best-fit nonlinear parameters"
            )
        self.keys["fit.amplitudes"] = self.schema.addField(
            "fit.amplitudes", type="ArrayD", size=self.model.getAmplitudeDim(),
            doc="best-fit linear amplitudes"
            )
        # If we're doing a warm start from a previous run with a different tag, check that things
        # are compatible, and add a subtask that matches what was used to run the previous.
        # This can be recursive.
        if self.config.previous is not None:
            if butler is None:
                raise lsst.pipe.base.TaskError("Cannot use previous outputs for warm start without butler;"
                                               " make sure you have have pipe_base #3085 if running from"
                                               " the command line")
            prevConfig = self.getPreviousConfig(butler)
            PrevTask = self.getPreviousTaskClass()
            self.previous = PrevTask(name="previous", config=prevConfig, parentTask=self, butler=butler)
            if (self.previous.model.getNonlinearDim() != self.model.getNonlinearDim()
                or self.previous.model.getAmplitudeDim() != self.model.getAmplitudeDim()
                or self.previous.model.getFixedDim() != self.model.getFixedDim()):
                raise lsst.pipe.base.TaskError("Cannot use previous catalog: model is incompatible")
            self.prevCatMapper = lsst.afw.table.SchemaMapper(self.previous.schema)
            self.prevCatMapper.addMinimalSchema(self.schema)
            assert self.prevCatMapper.getOutputSchema() == self.schema
        else:
            self.previous = None
        # Create the fitter subtask that does all the non-bookkeeping work
        self.makeSubtask("fitter", schema=self.schema, keys=self.keys, model=self.model, prior=self.prior,
                         previous=(self.previous.fitter if self.previous else None))

    def makeUnitSystem(self, record, position, magnitude):
        record.set(self.keys['sys.position'], position)
        record.set(self.keys['sys.magnitude'], magnitude)
        return modelfitLib.UnitSystem(position, magnitude)

    def getUnitSystem(self, record):
        return modelfitLib.UnitSystem(record.get(self.keys['sys.position']),
                                      record.get(self.keys['sys.magnitude']))

    def makeTable(self):
        """Return a ModelFitTable object based on the measurement schema and the fitter subtask's
        sample schema.
        """
        table = lsst.meas.modelfit.ModelFitTable.make(self.schema, self.fitter.makeSampleTable(),
                                                      self.fitter.interpreter)
        return table

    def adaptPrevious(self, prevCat):
        """Adapt a previous catalog to create a new output catalog."""
        outCat = modelfitLib.ModelFitCatalog(self.makeTable())
        outCat.extend(prevCat, mapper=self.prevCatMapper)
        for n, (prevRecord, outRecord) in enumerate(zip(prevCat, outCat)):
            if self.config.progressChunk > 0 and n % self.config.progressChunk == 0:
                self.log.info("Adapting objects %d-%d of %d (currently %3.2f%%)"
                              % (n+1, n+self.config.progressChunk, len(outCat), (100.0*n)/len(outCat)))
            self.fitter.adaptPrevious(prevRecord, outRecord)
        return outCat

    def run(self, dataRef):
        """Main driver for model fitting.

        Subclasses should not need to override this method, and instead should customize
        readInputs(), prepCatalog(), makeLikelihood(), and writeOutputs().
        """
        self.log.info("Reading inputs")
        inputs = self.readInputs(dataRef)
        if self.previous is not None:
            self.log.info("Using previous run with tag '%s' for warm start" % self.config.previous)
            outCat = self.adaptPrevious(inputs.prevCat)
        else:
            self.log.info("Preparing catalog")
            outCat = self.prepCatalog(inputs)
        if self.config.maxObjects is not None:
            outCat = outCat[:self.config.maxObjects]
        if not self.config.prepOnly:
            for n, outRecord in enumerate(outCat):
                if self.config.progressChunk > 0 and n % self.config.progressChunk == 0:
                    self.log.info("Fitting objects %d-%d of %d (currently %3.2f%%)"
                                  % (n+1, n+self.config.progressChunk, len(outCat), (100.0*n)/len(outCat)))
                likelihood = self.makeLikelihood(inputs, outRecord)
                try:
                    self.fitter.run(likelihood, outRecord)
                    nonlinear = outRecord[self.keys['fit.nonlinear']]
                    nonlinear[:] = self.fitter.interpreter.computeNonlinearMean(outRecord)
                    amplitudes = outRecord[self.keys['fit.amplitudes']]
                    amplitudes[:] = self.fitter.interpreter.computeAmplitudeMean(outRecord)
                except Exception as err:
                    if self.config.doRaise:
                        raise
                    else:
                        self.log.warn("Failure fitting object %d of %d with ID=%d:"
                                      % (n, len(outCat), outRecord.getId()))
                        self.log.warn(str(err))
                        continue
        self.log.info("Writing output catalog")
        self.writeOutputs(dataRef, outCat)
        return lsst.pipe.base.Struct(outCat=outCat, inputs=inputs)

    def getPreviousTaskClass(self):
        """Return the Task class used to construct the previous catalog, if applicable."""
        raise NotImplementedError()

    def getPreviousConfig(self, butler):
        """Fetch and return the config instance used to construct the previous catalog, if applicable."""
        raise NotImplementedError()

    def readInputs(self, dataRef):
        """Read task inputs using the butler.

        The form of the return value is subclass-dependent; it will simply be passed unmodified
        to prepCatalog() and makeLikelihood()
        """
        raise NotImplementedError()

    def prepCatalog(self, inputs):
        """Prepare and return the output catalog, doing everything but the actual fitting.

        After this step, each output record should be in a state such that makeLikelihood and
        fitter.run() may be called on it.
        """
        raise NotImplementedError()

    def makeLikelihood(self, inputs, record):
        """Create a Likelihood object for a single object.

        The inputs are passed unmodified from the readInputs() method, and the record is
        as produced by prepCatalog().
        """
        raise NotImplementedError()

    def writeOutputs(self, dataRef, outCat):
        """Write task outputs using the butler.
        """
        raise NotImplementedError()

    def writeConfig(self, butler, clobber=False):
        dataRef = butler.dataRef(self._getConfigName(), tag=self.config.tag)
        # We trick the base class implementation into including a data ID
        # by giving it a DataRef instead of a Butler.
        lsst.pipe.base.CmdLineTask.writeConfig(self, butler=dataRef, clobber=clobber)

    def writeSchemas(self, butler, clobber=False):
        dataRef = butler.dataRef(self._getConfigName(), tag=self.config.tag)
        # We trick the base class implementation into including a data ID
        # by giving it a DataRef instead of a Butler.
        lsst.pipe.base.CmdLineTask.writeSchemas(self, butler=dataRef, clobber=clobber)

    def writeMetadata(self, dataRef):
        # Reimplement CmdLineTask.writeMetadata to include tag in dataId
        dataRef.put(self.getFullMetadata(), self._getMetadataName(), tag=self.config.tag)
