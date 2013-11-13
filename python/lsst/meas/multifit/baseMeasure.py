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
        default="bulge+disk",
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
    fitPixelScale = lsst.pex.config.Field(
        dtype=float,
        default=0.2,
        doc="Pixel scale (arcseconds/pixel) for coordinate system used for model parameters"
    )
    fitFluxMag0 = lsst.pex.config.Field(
        dtype=float,
        default=30.0,
        doc="Flux at magnitude 0 used for to define the units of amplitude in models"
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

    def makeFitWcs(self, coord):
        return lsst.afw.image.makeLocalWcs(coord, self.fitPixelScale * lsst.afw.geom.arcseconds)

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
        self.prior = self.config.prior.apply(pixelScale=self.config.fitPixelScale*lsst.afw.geom.arcseconds,
                                             fluxMag0=self.config.fitFluxMag0)
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
        # This Calib determines the flux units we use for amplitude parameters; like the WCS,
        # we want this to be a global system so all Objects have the same units
        self.fitCalib = lsst.afw.image.Calib()
        self.fitCalib.setFluxMag0(self.config.fitFluxMag0)
        # Create the fitter subtask that does all the non-bookkeeping work
        self.makeSubtask("fitter", schema=self.schema, keys=self.keys, model=self.model, prior=self.prior)

    def makeTable(self):
        """Return a ModelFitTable object based on the measurement schema and the fitter subtask's
        sample schema.
        """
        table = lsst.meas.multifit.ModelFitTable.make(self.schema, self.fitter.makeTable())
        table.setInterpreter(self.fitter.interpreter)
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
                    self.log.info("Processing object %d/%d (%3.2f%%)"
                                  % (n, len(outCat), (100.0*n)/len(outCat)))
                likelihood = self.makeLikelihood(inputs, outRecord)
                try:
                    self.fitter.run(likelihood, outRecord)
                except Exception as err:
                    self.log.warn("Failure fitting object %d of %d with ID=%d"
                                  % (n, len(outCat), outRecord.getId()))
                    continue
        self.log.info("Writing output catalog")
        self.writeOutputs(dataRef, outCat)
        return outCat

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
