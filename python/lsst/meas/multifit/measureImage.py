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

from .samplers import BaseSamplerTask, AdaptiveImportanceSamplerTask
from .multifitLib import ProjectedLikelihood, ModelFitCatalog, ModelFitTable
from .fitRegion import setupFitRegion
from .models import *
from .priors import *

__all__ = ("BaseMeasureConfig", "MeasureImageConfig", "MeasureImageTask")

class BaseMeasureConfig(lsst.pex.config.Config):
    sampler = lsst.pex.config.ConfigurableField(
        target=AdaptiveImportanceSamplerTask,
        doc="Subtask that generates samples from the probability of a galaxy model given image data"
    )
    model = modelRegistry.makeField(
        default="bulge+disk",
        doc="Definition of the galaxy model to fit"
    )
    prior = priorRegistry.makeField(
        default="mixture",
        doc="Bayesian prior on galaxy parameters"
    )
    psf = lsst.pex.config.ConfigField(
        dtype=lsst.meas.extensions.multiShapelet.FitPsfConfig,
        doc="Config options for approximating the PSF using shapelets"
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
    prepOnly = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        doc="If True, only prepare the catalog (match, transfer fields, fit PSF)"
    )
    marginalizeAmplitudes = lsst.pex.config.Field(
        dtype=bool,
        default=True,
        doc="Marginalize over amplitudes numerically instead of sampling them?"
    )

    def setDefaults(self):
        self.psf.innerOrder = 4
        self.psf.outerOrder = 0

class MeasureImageConfig(BaseMeasureConfig):
    likelihood = lsst.pex.config.ConfigField(
        dtype=ProjectedLikelihood.ConfigClass,
        doc="Config for likelihood object that computes model probability at given parameters"
    )
    doWarmStart = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        doc=("If True, load a previous modelfits catalog and use its attached proposal distributions "
             "as the initial proposal distributions, instead of starting from a match of the source "
             "and reference catalogs.  NOTE: Also causes fit region footprints to be based on the previous "
             "modelfits footprints, instead of the original detection footprints")
    )

class MeasureImageTask(lsst.pipe.base.CmdLineTask):
    """Driver class for S13-specific galaxy modeling work

    Like ProcessImageTask, MeasureImageTask is intended to be used as a base
    class with CCD and coadd derived-class specializations.  It is run
    after processCcd.py (or processEImage.py), and generates a single output
    catalog with the mapper name 'modelfits'.
    """

    ConfigClass = MeasureImageConfig
    dataPrefix = ""

    def __init__(self, **kwds):
        BaseMeasureTask.__init__(self, **kwds)
        self.schemaMapper = lsst.afw.table.SchemaMapper(lsst.afw.table.SourceTable.makeMinimalSchema())
        self.schemaMapper.addMinimalSchema(lsst.meas.multifit.ModelFitTable.makeMinimalSchema())
        self.schema = self.schemaMapper.getOutputSchema()
        self.sampleSchema = lsst.afw.table.Schema()
        self.objectiveFactory = self.config.objectiveFactory(
        self.makeSubtask("sampler", sampleSchema=)
        self.fitPsf = self.config.psf.makeControl().makeAlgorithm(self.schema)
        self.model = self.config.model.apply()
        self.prior = self.config.prior.apply(self.config.fitPixelScale, self.config.fitFluxMag0)
        self.calib = lsst.afw.image.Calib()
        self.calib.setFluxMag0(self.config.fitFluxMag0)
        self.keys = {}
        self.keys["ref.center"] = self.schema.addField(
            "ref.center", type="PointD",
            doc="position in image coordinates from reference catalog"
            )
        self.keys["ref.parameters"] = self.schema.addField(
            "ref.parameters", type="ArrayD", size=self.model.getParameterDim(),
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
        self.keys["mean.parameters"] = self.schema.addField(
            "mean.parameters", type="ArrayD", size=self.model.getParameterDim(),
            doc="posterior mean nonlinear parameters"
            )

    def readInputs(self, dataRef):
        """Return a lsst.pipe.base.Struct containing:
          - exposure ----- lsst.afw.image.ExposureF to fit
          - srcCat ------- lsst.afw.table.SourceCatalog with initial measurements
          - refCat ------- lsst.afw.table.SimpleCatalog with truth values
        """
        return lsst.pipe.base.Struct(
            exposure = dataRef.get(self.dataPrefix + "calexp", immediate=True),
            srcCat = dataRef.get(self.dataPrefix + "src", immediate=True),
            refCat = dataRef.get("refcat", immediate=True)
            )

    def writeOutputs(self, dataRef, outCat):
        """Write task outputs using the butler.
        """
        self.log.info("Writing output catalog")
        dataRef.put(outCat, self.dataPrefix + "modelfits")

    @lsst.pipe.base.timeMethod
    def processObject(self, exposure, record, doWarmStart=False):
        """Process a single object.

        @param[in] exposure     lsst.afw.image.ExposureF to fit
        @param[in,out] record   ModelFitRecord to fill, as prepared by prepCatalog
        @param[in] doWarmStart  If True, attempt to load the initial proposal distribution
                                from a SampleSet already attached to the given record.

        @return a Struct containing various intermediate objects and results:
          - likelihood: the Likelihood object used to evaluate model probabilities
          - sampler: the Sampler object used to draw samples
          - psf: a shapelet.MultiShapeletFunction representation of the PSF
          - record: the output record (identical to the record argument, which is modified in-place)
        """
        if record.getSchema() != self.schema:
            raise TaskError("Record schema does not match expected schema; probably a result of using "
                            "doWarmStart with a catalog from an older version of the code.")
        psfModel = lsst.meas.extensions.multiShapelet.FitPsfModel(self.config.psf.makeControl(), record)
        psf = psfModel.asMultiShapelet()
        center = record.getPointD(self.keys["ref.center"])
        if not doWarmStart:
            sampler = self.sampler.setup(exposure=exposure, center=center, prior=self.prior,
                                         ellipse=record.getMomentsD(self.keys["ref.ellipse"]))
        else:
            samples = record.getSamples()
            if samples is None:
                raise TaskError("No prior samples found; cannot proceed with warm start")
            sampler = self.sampler.reset(samples=samples, center=center, prior=self.prior)
        likelihood = SingleEpochLikelihood(
            self.config.likelihood.makeControl(), self.basis, psf,
            exposure.getMaskedImage(), record.getFootprint()
        )
        samples = sampler.run(likelihood)
        samples.applyPrior(self.prior)
        record.setSamples(samples)
        self.fillDerivedFields(record)
        return lsst.pipe.base.Struct(likelihood=likelihood, sampler=sampler, psf=psf, record=record)

    def prepCatalog(self, exposure, srcCat, refCat, where=None):
        """Create a ModelFitCatalog with initial parameters and fitting regions
        to be used later by fillCatalog.

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
        self.log.info("Setting up fitting and transferring source/reference fields")
        outCat = ModelFitCatalog(self.schema)
        if not exposure.hasPsf():
            raise lsst.pipe.base.TaskError("Exposure has no PSF")
        matches = lsst.afw.table.matchRaDec(refCat, srcCat, 1.0*lsst.afw.geom.arcseconds)
        keyA = refCat.getSchema().find("ellipse.a").key
        keyB = refCat.getSchema().find("ellipse.b").key
        keyTheta = refCat.getSchema().find("ellipse.theta").key
        keyMag = refCat.getSchema().find("mag.%s" % exposure.getFilter().getName()).key
        keySIndex = refCat.getSchema().find("sindex").key

        for match in matches:
            srcRecord = match.second
            refRecord = match.first
            if where is not None and not where(src=srcRecord, ref=refRecord):
                continue
            if srcRecord.getShapeFlag():
                continue
            outRecord = outCat.addNew()
            outRecord.assign(srcRecord, self.schemaMapper)
            wcs = lsst.afw.image.makeLocalWcs(
                refRecord.getCoord(), self.config.fitPixelScale * lsst.afw.geom.arcseconds
                )

            ellipse1 = lsst.afw.geom.ellipses.Ellipse(
                lsst.afw.geom.ellipses.Axes(
                    (refRecord.getD(keyA)*lsst.afw.geom.arcseconds).asDegrees(),
                    (refRecord.getD(keyB)*lsst.afw.geom.arcseconds).asDegrees(),
                    (refRecord.getD(keyTheta)*lsst.afw.geom.degrees).asRadians() - 0.5*numpy.pi
                    ),
                refRecord.getCoord().getPosition(lsst.afw.geom.degrees)
                )
            transform = wcs.linearizeSkyToPixel(refRecord.getCoord())
            ellipse2 = ellipse1.transform(transform)
            ellipses = lsst.meas.multifit.Model.EllipseVector([ellipse2])
            self.model.readEllipses(
                ellipses,
                outRecord[self.keys["ref.parameters"]],
                outRecord[self.keys["ref.fixed"]]
                )
            # this flux->amplitudes conversion assumes ref catalog is single-component, and that the
            # first component of the model is what that corresponds to; we'll need to generalize
            # this eventually
            flux = self.calib.getFlux(refRecord.get(keyMag))
            outRecord[self.keys["ref.amplitudes"]][:] = 0.0
            outRecord.setD(self.keys["ref.amplitudes"][0], flux)
            outRecord.setD(self.keys["snr"], srcRecord.getApFlux() / srcRecord.getApFluxErr())
            outRecord.setPointD(self.keys["ref.center"], exposure.getWcs().skyToPixel(refRecord.getCoord()))
            outRecord.setFootprint(setupFitRegion(self.config.fitRegion, exposure, srcRecord))
            self.fitPsf.fit(outRecord, exposure.getPsf(), outRecord.get(self.keys["ref.center"]))
        outCat.sort()
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
        if self.config.doWarmStart:
            self.log.info("Using warm start; loading previous catalog and updating footprints and PSFs")
            outCat = self.dataRef.get(self.dataPrefix + "src", immediate=True)
            for outRecord in outCat:
                outRecord.setFootprint(setupFitRegion(self.config.fitRegion, inputs.exposure, outRecord))
                psfModel = self.fitPsf.fit(outRecord, inputs.exposure.getPsf(),
                                           outRecord.get(self.keys["ref.center"]))
        else:
            outCat = self.prepCatalog(inputs.exposure, srcCat=inputs.srcCat, refCat=inputs.refCat)
        if not self.config.prepOnly:
            self.fillCatalog(inputs.exposure, outCat)
        self.writeOutputs(dataRef, outCat)
        return outCat

    def getSchemaCatalogs(self):
        """Return a dict of empty catalogs for each catalog dataset produced by this task."""
        return {self.dataPrefix + "modelfits": ModelFitCatalog(self.schema)}
