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
import sys
import traceback

import numpy

from lsst.pipe.base import Struct, TaskError
import lsst.pex.config
import lsst.afw.geom as afwGeom
from lsst.meas.extensions.multiShapelet import FitPsfAlgorithm
from .multifitLib import VectorEpochFootprint, EpochFootprint, MultiEpochObjective, ModelFitCatalog, \
    ModelFitTable
from .measureImage import BaseMeasureConfig, BaseMeasureTask
from .samplers import *

try:
    from lsst.meas.mosaic import applyMosaicResults
except ImportError:
    applyMosaicResults = None

__all__ = ("MeasureMultiConfig", "MeasureMultiTask")

class MeasureMultiConfig(BaseMeasureConfig):
    coaddName = lsst.pex.config.Field(
        doc = "coadd name: typically one of deep or goodSeeing",
        dtype = str,
        default = "deep",
    )
    objective = lsst.pex.config.ConfigField(
        dtype=MultiEpochObjective.ConfigClass,
        doc="Config for objective object that computes model probability at given parameters"
    )
    minPixels = lsst.pex.config.Field(
        doc = "minimum number of pixels in a calexp footprint to use that calexp for a given galaxy",
        dtype = int,
        default = 5,
    )
    doApplyUberCal = lsst.pex.config.Field(
        dtype = bool,
        doc = "Apply meas_mosaic ubercal results to input calexps?",
        default = True
    )

    def setDefaults(self):
        self.sampler.iterations.clear()
        self.sampler.iterations[0] = ImportanceSamplerConfig(nUpdateSteps=0, nSamples=200,
                                                             maxRepeat=0, targetPerplexity=1.0)

class MeasureMultiTask(BaseMeasureTask):
    """Variant of MeasureImageTask for running multifit on the calexp that make up a coadd.

    The tasks are so different in implementation that no code is shared (yet).
    """
    ConfigClass = MeasureMultiConfig
    _DefaultName = "measureMulti"

    def __init__(self, **kwds):
        BaseMeasureTask.__init__(self, **kwds)
        self.dataPrefix = self.config.coaddName + "Coadd_"
        self.makeSubtask("sampler")
        self.schema = None # set later when we have a coadd catalog to copy
        self.basis = self.config.model.apply()
        self.prior = None # set later when we can look up coadd pixel scale
        self.psfControl = self.config.psf.makeControl()
        self.keys = {}
        self.readInputExposure = None  # placeholder for closure to be installed by readInputs

    def run(self, dataRef):
        """Process the catalog associated with the given coadd

        For each source in the catalog find those coadd input images that fully overlap the footprint,
        and fit the source shape. Combine the fits from all input images to form a multifit catalog.

        @param[in] dataRef: coadd data reference
        """
        inputs = self.readInputs(dataRef)
        outCat = self.prepCatalog(coaddCat=inputs.coaddCat)
        if not self.config.prepOnly:
            numObjects = len(outCat)
            for i, record in enumerate(outCat):
                self.log.info("Processing object %s of %s" % (i + 1, numObjects))
                try:
                    self.processObject(record=record, coadd=inputs.coadd, coaddInputCat=inputs.coaddInputCat)
                except Exception, e:
                    self.log.warn("processObject failed: %s\n" % (e,))
                    traceback.print_exc(file=sys.stderr)
        self.writeOutputs(dataRef, outCat)
        return outCat

    def readInputs(self, dataRef):
        """Return inputs, and attach to the task a closure method (readInputExposure)
        that loads a calexp subimage given an ExposureRecord and a bounding box.

        @param[in] dataRef: data reference for coadd

        @return an lsst.pipe.base.Struct containing:
          - coadd: coadd patch (lsst.afw.image.ExposureF); the images that make up this coadd are fit
          - coaddInputCat: catalog of ExposureRecords corresponding to calexps that went into the coadd
          - coaddCat: catalog with model fits based on the coadd (lsst.meas.multifit.ModelFitCatalog)
        """
        coadd = dataRef.get(self.dataPrefix[0:-1], immediate=True)
        coaddInputCat = coadd.getInfo().getCoaddInputs().ccds
        coaddInputSchema = coaddInputCat.getSchema()
        visitKey = coaddInputSchema.find("visit").key
        ccdKey = coaddInputSchema.find("ccd").key
        butler = dataRef.getButler()
        def readInputExposure(record, bbox):
            """Given an ExposureRecord and bounding box, load the appropriate subimage."""
            dataId = butler.mapper.getDataId(visit=record.get(visitKey), ccdId=record.get(ccdKey))
            dataRef = butler.dataRef("calexp", **dataId)
            exposure = dataRef.get(
                "calexp_sub",
                bbox=bbox,
                origin="PARENT",
                immediate=True
            )
            if self.config.doApplyUberCal:
                if not applyMosaicResults:
                    raise RuntimeError(
                        "Cannot use improved calibrations for %s because meas_mosaic could not be imported."
                        % dataRef.dataId
                        )
                applyMosaicResults(dataRef, calexp=exposure, bbox=bbox)
            return exposure
        self.readInputExposure = readInputExposure
        return lsst.pipe.base.Struct(
            coadd=coadd,
            coaddInputCat=coaddInputCat,
            coaddCat=dataRef.get(self.dataPrefix + "modelfits", immediate=True),
        )

    def prepCatalog(self, coaddCat):
        """Create an output ModelFitCatalog that copies appropriate fields from the coadd ModelFitCatalog

        @param[in] exposure     Exposure object that will be fit.
        @param[in] coaddCat     SourceCatalog containing MeasureCoaddTask measurements
        @return a ModelFit catalog with one entry per coaddCat
        """
        if self.schema is None:
            self.setSchema(coaddCat)
        self.log.info("Copying data from coadd catalog to new catalog")
        outCat = ModelFitCatalog(self.schema)
        for coaddRecord in coaddCat:
            record = outCat.addNew()
            record.assign(coaddRecord, self.schemaMapper)
        outCat.sort()
        return outCat

    @lsst.pipe.base.timeMethod
    def processObject(self, record, coadd, coaddInputCat, doWarmStart=True):
        """Process a single object.

        @param[in,out] record     multi-fit ModelFitRecord
        @param[in] coadd          Coadd exposure
        @param[in] coaddInputCat    ExposureCatalog describing exposures that may be included in the fit
        @param[in] doWarmStart      If True (default), use the proposal distribution attached to the
                                    coadd-fit SampleSet as the initial proposal.

        @return a Struct containing various intermediate objects:
          - objective   the Objective object used to evaluate likelihoods
          - sampler     the Sampler object used to draw samples
          - record      the output record (identical to the record argument, which is modified in-place)
        """
        if self.prior is None:
            self.prior = self.config.prior.apply(pixelScale=coadd.getWcs().pixelScale())
        objective = self.makeObjective(
            record=record,
            coadd=coadd,
            coaddInputCat=coaddInputCat,
        )
        refCenter = record.getPointD(self.keys["ref.center"])
        if doWarmStart:
            samples = record.getSamples()
            if samples is None:
                raise TaskError("No prior samples found; cannot proceed with warm start")
            sampler = self.sampler.reset(samples=samples, center=refCenter, prior=self.prior)
        else:
            sampler = self.sampler.setup(
                exposure=coadd,
                center=refCenter,
                ellipse=record.getMomentsD(self.keys["ref.ellipse"]),
                )
        samples = sampler.run(objective)
        samples.applyPrior(self.prior)
        record.setSamples(samples)
        self.fillDerivedFields(record)
        return lsst.pipe.base.Struct(objective=objective, sampler=sampler, record=record)

    @lsst.pipe.base.timeMethod
    def makeObjective(self, record, coadd, coaddInputCat):
        """Construct a MultiEpochObjective from the calexp footprints

        @param[in] record       multi-fit ModelFitRecord
        @param[in] coadd        Coadd exposure
        @param[in] coaddInputCat  ExposureCatalog describing exposures that may be included in the fit
        """
        if self.keys is None:
            self.keys = {"source.center": record.getSchema().find("source.center").key}
        coaddFootprint = record.getFootprint()
        sourceCoaddPos = record.getPointD(self.keys["source.center"])

        coaddWcs = coadd.getWcs()
        sourceSkyPos = coaddWcs.pixelToSky(sourceCoaddPos)

        # process each calexp that partially overlaps this footprint
        epochFootprintList = VectorEpochFootprint()

        for exposureRecord in coaddInputCat:
            calexpFootprint = coaddFootprint.transform(coaddWcs, exposureRecord.getWcs(),
                                                       exposureRecord.getBBox())
            calexpFootprint.clipTo(exposureRecord.getBBox()) # Footprint.transform does not clip

            # due to ticket #2979 this test is invalid -- getNpix includes clipped pixels!
            # so also test that the footprint's bbox is not empty
            if calexpFootprint.getNpix() < self.config.minPixels:
                # no overlapping pixels, so skip this calexp
                continue
            calexpFootprintBBox = calexpFootprint.getBBox()
            if calexpFootprintBBox.isEmpty(): # temporary hack due to ticket #2979
                continue

            calexp = self.readInputExposure(record=exposureRecord, bbox=calexpFootprintBBox)

            sourceCalexpPos = calexp.getWcs().skyToPixel(sourceSkyPos)

            psfModel = FitPsfAlgorithm.apply(self.psfControl, calexp.getPsf(), sourceCalexpPos)
            psf = psfModel.asMultiShapelet()

            epochFootprint = EpochFootprint(calexpFootprint, calexp, psf)
            epochFootprintList.append(epochFootprint)

        return MultiEpochObjective(
            self.config.objective.makeControl(),
            self.basis,
            coaddWcs,
            sourceSkyPos,
            epochFootprintList,
        )

    def setSchema(self, coaddCat):
        """Construct self.schema; call once as soon as you have your first coadd catalog

        @param[in] coaddCat  ModelFitCatalog from coadd
        @raise RuntimeError if self.schema is not None
        """
        if self.schema is not None:
            raise RuntimeError("self.schema already set")

        self.log.info("Setting the schema")
        self.schemaMapper = lsst.afw.table.SchemaMapper(coaddCat.getSchema())
        self.schemaMapper.addMinimalSchema(ModelFitTable.makeMinimalSchema())

        def mapKey(keyName):
            """Map one key from coaddCat to self.schemaMapper
            @param[in] keyName      name of key to map
            """
            inputKey = coaddCat.getSchema().find(keyName).getKey()
            self.keys[keyName] = self.schemaMapper.addMapping(inputKey)

        for name in ("source.ellipse", "source.center", "snr", "ref.ellipse", "ref.center", "ref.sindex"):
            mapKey(name)

        self.schema = self.schemaMapper.getOutputSchema()
        self.addDerivedFields()

    def writeOutputs(self, dataRef, outCat):
        """Write task outputs using the butler.
        """
        self.log.info("Writing output catalog")
        dataRef.put(outCat, self.dataPrefix + "multiModelfits")

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser
        """
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "deepCoadd",
            help="coadd data ID, e.g. --id tract=1 patch=2,2 filter=g")
        return parser

    def _getConfigName(self):
        """Return the name of the config dataset
        """
        return "%s_measureMulti_config" % (self.config.coaddName,)

    def _getMetadataName(self):
        """Return the name of the metadata dataset
        """
        return "%s_measureMulti_metadata" % (self.config.coaddName,)
