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
import lsst.afw.geom.ellipses
import lsst.meas.extensions.multiShapelet

from .baseMeasure import *
from .multifitLib import SingleEpochObjective, ModelFitCatalog, ModelFitTable

__all__ = ("MeasureImageConfig", "MeasureImageTask")

class MeasureImageConfig(InitialMeasureConfig):
    objective = lsst.pex.config.ConfigField(
        dtype=SingleEpochObjective.ConfigClass,
        doc="Config for objective object that computes model probability at given parameters"
    )
    useRefCat = lsst.pex.config.Field(
        dtype=bool,
        default=True,
        doc="Whether to use the reference catalog to identify objects to fit"
    )

class MeasureImageTask(InitialMeasureTask):
    """Driver class for S13-specific galaxy modeling work

    Like ProcessImageTask, MeasureImageTask is intended to be used as a base
    class with CCD and coadd derived-class specializations.  It is run
    after processCcd.py (or processEImage.py), and generates a single output
    catalog with the mapper name 'modelfits'.
    """

    ConfigClass = MeasureImageConfig
    dataPrefix = ""

    def __init__(self, **kwds):
        InitialMeasureTask.__init__(self, **kwds)
        self.makeSubtask("sampler")
        self.basis = self.config.model.apply()
        self.fitPsf = None
        self.prior = None

    def readInputs(self, dataRef):
        """Return a lsst.pipe.base.Struct containing:
          - exposure ----- lsst.afw.image.ExposureF to fit
          - srcCat ------- lsst.afw.table.SourceCatalog with initial measurements
          - refCat ------- lsst.afw.table.SimpleCatalog with truth values
        """
        try:
            refCat = dataRef.get("refcat", immediate=True)
        except:
            refCat = None
        if refCat is None and self.config.useRefCat:
            raise TaskError("useRefCat true but no reference catalog found")
        return lsst.pipe.base.Struct(
            exposure = dataRef.get(self.dataPrefix + "calexp", immediate=True),
            srcCat = dataRef.get(self.dataPrefix + "src", immediate=True),
            refCat = refCat
            )

    def writeOutputs(self, dataRef, outCat):
        """Write task outputs using the butler.
        """
        self.log.info("Writing output catalog")
        outCat.sort()  # want to sort by ID before saving, so when we load it's contiguous
        dataRef.put(outCat, self.dataPrefix + "modelfits")

    @lsst.pipe.base.timeMethod
    def processObject(self, exposure, record):
        """Process a single object.

        @param[in] exposure    lsst.afw.image.ExposureF to fit
        @param[in,out] record  ModelFitRecord to fill, as prepared by prepCatalog

        @return a Struct containing various intermediate objects and results:
          - objective: the Objective object used to evaluate likelihoods
          - sampler: the Sampler object used to draw samples
          - psf: a shapelet.MultiShapeletFunction representation of the PSF
          - record: the output record (identical to the record argument, which is modified in-place)
        """
        if self.prior is None:
            self.prior = self.config.prior.apply(pixelScale=exposure.getWcs().pixelScale())
        if not self.schema: # we may not have run prepCatalog, in which case we need some keys
            self.schema = record.getSchema()
            self.keys["ref.center"] = self.schema.find("ref.center").key
            self.keys["ref.ellipse"] = self.schema.find("ref.ellipse").key
            self.keys["source.center"] = self.schema.find("source.center").key
            self.keys["source.ellipse"] = self.schema.find("source.ellipse").key
            self.keys["mean.center"] = self.schema.find("mean.center").key
            self.keys["mean.ellipse"] = self.schema.find("mean.ellipse").key
            self.keys["median.center"] = self.schema.find("median.center").key
            self.keys["median.ellipse"] = self.schema.find("median.ellipse").key
        psfModel = lsst.meas.extensions.multiShapelet.FitPsfModel(self.config.psf.makeControl(), record)
        psf = psfModel.asMultiShapelet()
        sampler = self.sampler.setup(exposure=exposure, center=record.getPointD(self.keys["ref.center"]),
                                     ellipse=record.getMomentsD(self.keys["ref.ellipse"]), prior=self.prior)
        objective = SingleEpochObjective(
            self.config.objective.makeControl(), self.basis, psf,
            exposure.getMaskedImage(), record.getFootprint()
        )
        samples = sampler.run(objective)
        samples.applyPrior(self.prior)
        record.setSamples(samples)
        self.fillMeasuredFields(record)
        return lsst.pipe.base.Struct(objective=objective, sampler=sampler, psf=psf, record=record)

    def setupInitialSchema(self, exposure, srcCat, refCat):
        """Attach Schema and SchemaMapper attributes to self (.schema and self.srcMapper, respectively)
        for use with prepInitialCatalog().

        MeasureImageTask simply uses the implementation in its base class and adds fields for the
        shapelet PSF fit (while initializing the self.fitPsf attribute).
        """
        InitialMeasureTask.setupInitialSchema(self, exposure=exposure, srcCat=srcCat, refCat=refCat)
        self.fitPsf = self.config.psf.makeControl().makeAlgorithm(self.schema)

    def prepInitialCatalog(self, exposure, srcCat, refCat, where=None):
        """Create a ModelFitCatalog with initial parameters and fitting regions
        to be used later by fillCatalog.

        MeasureImageTask simply uses the implementation in its base class and adds the
        shapelet PSF fit.

        @param[in]   exposure         Exposure object that will be fit.
        @param[in]   srcCat           SourceCatalog containing processCcd measurements
        @param[in]   refCat           Simulation reference catalog used to add truth
                                      values for comparison to each record, and to
                                      set which objects should be fit (i.e. reject stars).
        @param[in]   where            Callable with signature f(src, ref=None) that takes
                                      a SourceRecord and optional SimpleRecord and returns
                                      True if a ModelFitRecord should be created for that
                                      object.  Ignored if None.
        """
        outCat = InitialMeasureTask.prepCatalog(self, exposure=exposure, srcCat=srcCat, refCat=refCat,
                                                where=where)
        self.log.info("Fitting shapelet PSF approximations")
        if not exposure.hasPsf():
            raise lsst.pipe.base.TaskError("Exposure has no PSF")
        for outRecord in outCat:
            self.fitPsf.fit(outRecord, exposure.getPsf(), srcRecord.getCentroid())
        return outCat

    def fillCatalog(self, exposure, outCat):
        """For each empty ModelFitRecord in catalog, call processObject()
        """
        for n, record in enumerate(outCat):
            if self.config.progressChunk > 0 and n % self.config.progressChunk == 0:
                self.log.info("Processing object %d/%d (%3.2f%%)" % (n, len(outCat), (100.0*n)/len(outCat)))
            self.processObject(exposure=exposure, record=record)

    def run(self, dataRef):
        """Process the exposure/catalog associated with the given dataRef, and write the
        output catalog using the butler.
        """
        inputs = self.readInputs(dataRef)
        outCat = self.prepCatalog(inputs.exposure, srcCat=inputs.srcCat, refCat=inputs.refCat)
        if not self.config.prepOnly:
            self.fillCatalog(inputs.exposure, outCat)
        self.writeOutputs(dataRef, outCat)
        return outCat

    def getSchemaCatalogs(self):
        """Return a dict of empty catalogs for each catalog dataset produced by this task."""
        return {self.dataPrefix + "modelfits": ModelFitCatalog(self.schema)}
