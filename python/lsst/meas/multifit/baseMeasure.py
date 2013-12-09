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
import lsst.meas.extensions.multiShapelet

from . import multifitLib
from .models import modelRegistry
from .priors import priorRegistry
from .samplers import AdaptiveImportanceSamplerTask
from .fitRegion import setupFitRegion

__all__ = ("BaseMeasureConfig", "BaseMeasureTask")

class BaseMeasureConfig(lsst.pex.config.Config):
    fitter = lsst.pex.config.ConfigurableField(
        target=AdaptiveImportanceSamplerTask,
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
        dtype=lsst.meas.extensions.multiShapelet.FitPsfConfig,
        doc="Config options for approximating the PSF using shapelets"
    )
    prepOnly = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        doc="If True, only prepare the catalog (match, transfer fields, fit PSF)"
    )
    doRaise = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        doc="If True, raise exceptions when individual objects fail instead of warning."
    )

class BaseMeasureTask(lsst.pipe.base.CmdLineTask):
    """An intermediate base class for top-level model-fitting tasks.

    Subclasses (for different kinds of input data) must implement three methods:
    readInputs(), prepCatalog(), makeLikelihood(), and writeOutputs().

    Different fitting/sampling algorithms can be plugged in via the fitter subtask.
    """

    ConfigClass = BaseMeasureConfig

    def __init__(self, **kwds):
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
        self.schema = multifitLib.ModelFitTable.makeMinimalSchema()
        self.keys = {}
        self.keys["ref.center"] = self.schema.addField(
            "ref.center", type="PointD",
            doc="position in image coordinates from reference catalog"
            )
        self.keys["ref.nonlinear"] = self.schema.addField(
            "ref.nonlinear", type="ArrayD", size=self.model.getNonlinearDim(),
            doc="nonlinear parameters from reference catalog"
            )
        self.keys["ref.amplitudes"] = self.schema.addField(
            "ref.amplitudes", type="ArrayD", size=self.model.getAmplitudeDim(),
            doc="linear amplitudes from reference catalog"
            )
        self.keys["ref.fixed"] = self.schema.addField(
            "ref.fixed", type="ArrayD", size=self.model.getFixedDim(),
            doc="fixed nonlinear parameters from reference catalog"
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
        # Create the fitter subtask that does all the non-bookkeeping work
        self.makeSubtask("fitter", schema=self.schema, keys=self.keys, model=self.model, prior=self.prior)

    def makeUnitSystem(self, record, position, magnitude):
        record.set(self.keys['sys.position'], position)
        record.set(self.keys['sys.magnitude'], magnitude)
        return multifitLib.UnitSystem(position, magnitude)

    def getUnitSystem(self, record):
        return multifitLib.UnitSystem(record.get(self.keys['sys.position']),
                                      record.get(self.keys['sys.magnitude']))

    def makeTable(self):
        """Return a ModelFitTable object based on the measurement schema and the fitter subtask's
        sample schema.
        """
        table = lsst.meas.multifit.ModelFitTable.make(self.schema, self.fitter.makeTable())
        return table

    def run(self, dataRef):
        """Main driver for model fitting.

        Subclasses should not need to override this method, and instead should customize
        readInputs(), prepCatalog(), makeLikelihood(), and writeOutputs().
        """
        self.log.info("Reading inputs")
        inputs = self.readInputs(dataRef)
        self.log.info("Preparing catalog")
        outCat = self.prepCatalog(inputs)
        if not self.config.prepOnly:
            for n, outRecord in enumerate(outCat):
                if self.config.progressChunk > 0 and n % self.config.progressChunk == 0:
                    self.log.info("Processing objects %d-%d of %d (currently %3.2f%%)"
                                  % (n+1, n+self.config.progressChunk, len(outCat), (100.0*n)/len(outCat)))
                likelihood = self.makeLikelihood(inputs, outRecord)
                try:
                    self.fitter.run(likelihood, outRecord)
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
        is as produced by prepCatalog().
        """
        raise NotImplementedError()

    def writeOutputs(self, dataRef, outCat):
        """Write task outputs using the butler.
        """
        raise NotImplementedError()
