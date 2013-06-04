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

from .samplers import BaseSamplerTask, NaiveGridSamplerTask
from .multifitLib import SingleEpochObjective, ModelFitCatalog, ModelFitTable
from .fitRegion import setupFitRegion
from .models import *
from .priors import *

__all__ = ("MeasureImageConfig", "MeasureImageTask")

class MeasureImageConfig(lsst.pex.config.Config):
    sampler = lsst.pex.config.ConfigurableField(
        target=NaiveGridSamplerTask,
        doc="Subtask that generates samples from the probability of a galaxy model given image data"
    )
    objective = lsst.pex.config.ConfigField(
        dtype=SingleEpochObjective.ConfigClass,
        doc="Config for objective object that computes model probability at given parameters"
    )
    model = modelRegistry.makeField(
        default="bulge+disk",
        doc="Definition of the galaxy model to fit"
    )
    prior = priorRegistry.makeField(
        default="single-component",
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
    useRefCat = lsst.pex.config.Field(
        dtype=bool,
        default=True,
        doc="Whether to use the reference catalog to identify objects to fit"
    )
    prepOnly = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        doc="If True, only prepare the catalog (match, transfer fields, fit PSF)"
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
        lsst.pipe.base.CmdLineTask.__init__(self, **kwds)
        self.makeSubtask("sampler")
        self.schemaMapper = lsst.afw.table.SchemaMapper(lsst.afw.table.SourceTable.makeMinimalSchema())
        self.schemaMapper.addMinimalSchema(lsst.meas.multifit.ModelFitTable.makeMinimalSchema())
        self.schema = self.schemaMapper.getOutputSchema()
        self.fitPsf = self.config.psf.makeControl().makeAlgorithm(self.schema)
        self.basis = self.config.model.apply()
        self.prior = self.config.prior.apply()
        self.snrKey = self.schema.addField("snr", type=float,
                                           doc="signal to noise ratio from source apFlux/apFluxErr")
        self.meanEllipseKey = self.schema.addField("mean.ellipse", type="MomentsD",
                                                   doc="Posterior mean ellipse")
        self.meanCenterKey = self.schema.addField("mean.center", type="PointD",
                                                  doc="Posterior mean center position")
        self.meanFluxKey = self.schema.addField("mean.flux", type=float, doc="Posterior mean flux")
        self.meanFractionKey = self.schema.addField("mean.fraction", type="ArrayD", size=self.basis.getSize(),
                                                    doc="Posterior mean flux fraction in each component")
        if self.config.useRefCat:
            self.refEllipseKey = self.schema.addField("ref.ellipse", type="MomentsD",
                                                      doc="ellipse from reference catalog")
            self.refFluxKey = self.schema.addField("ref.flux", type=float, doc="flux from reference catalog")
            self.refSIndexKey = self.schema.addField("ref.sindex", type=float,
                                                     doc="Sersic index from reference catalog")

    def readInputs(self, dataRef):
        """Return a lsst.pipe.base.Struct containing:
          - exposure ----- lsst.afw.image.ExposureF to fit
          - srcCat ------- lsst.afw.table.SourceCatalog with initial measurements
          - refCat ------- lsst.afw.table.SimpleCatalog with truth values (may be
                           None if !config.useRefCat)
        """
        try:
            refCat = dataRef.get("refcat", immediate=True)
        except:
            refCat = None
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

        @return the Objective object used to do the fitting
        """
        psfModel = lsst.meas.extensions.multiShapelet.FitPsfModel(self.config.psf.makeControl(), record)
        psf = psfModel.asMultiShapelet()
        sampler = self.sampler.setup(exposure=exposure, record=record)
        objective = SingleEpochObjective(
            self.config.objective.makeControl(), self.basis, psf,
            exposure.getMaskedImage(), record.getFootprint()
        )
        samples = sampler.run(objective)
        samples.applyPrior(self.prior)
        record.setSamples(samples)
        ellipse = samples.interpret(samples.computeMean(), record.getCentroid())
        record.set(self.meanEllipseKey, lsst.afw.geom.ellipses.Quadrupole(ellipse.getCore()))
        record.set(self.meanCenterKey, ellipse.getCenter())
        # TODO: fill fraction, flux keys, also error fields
        return lsst.pipe.base.Struct(objective=objective, sampler=sampler, record=record)

    def prepCatalog(self, exposure, srcCat, refCat=None, where=None):
        """Create a ModelFitCatalog with initial parameters and fitting regions
        to be used later by fillCatalog.

        @param[in]   exposure         Exposure object that will be fit; only used to
                                      extract the Psf, Wcs, Calib, and bounding box.
        @param[in]   srcCat           SourceCatalog containing processCcd measurements
        @param[in]   refCat           Simulation reference catalog used to add truth
                                      values for comparison to each record, and to
                                      set which objects should be fit (i.e. reject stars).
                                      Ignored if !config.useRefCat.
        @param[in]   where            Callable with signature f(src, ref=None) that takes
                                      a SourceRecord and optional SimpleRecord and returns
                                      True if a ModelFitRecord should be created for that
                                      object.  Ignored if None.
        """
        self.log.info("Setting up fitting and transferring source/reference fields")
        outCat = ModelFitCatalog(self.schema)
        if not exposure.hasPsf():
            raise lsst.pipe.base.TaskError("Exposure has no PSF")

        def prepRecord(outRecord, srcRecord):
            psfModel = self.fitPsf.apply(outRecord, exposure.getPsf(), srcRecord.getCentroid())
            outRecord.setFootprint(setupFitRegion(self.config.fitRegion, exposure, srcRecord))
            outRecord.setD(self.snrKey, srcRecord.getApFlux() / srcRecord.getApFluxErr())
            outRecord.setCentroid(srcRecord.getCentroid())
            outRecord.setShape(srcRecord.getShape())
            outRecord.assign(srcRecord, self.schemaMapper)

        if self.config.useRefCat:
            matches = lsst.afw.table.matchRaDec(refCat, srcCat, 1.0*lsst.afw.geom.arcseconds)
            keyA = refCat.getSchema().find("ellipse.a").key
            keyB = refCat.getSchema().find("ellipse.b").key
            keyTheta = refCat.getSchema().find("ellipse.theta").key
            keyMag = refCat.getSchema().find("mag.%s" % exposure.getFilter().getName()).key
            keySIndex = refCat.getSchema().find("sindex").key
            wcs = exposure.getWcs()
            calib = exposure.getCalib()
            for match in matches:
                if where is not None and not where(src=match.second, ref=match.first):
                    continue
                ellipse1 = lsst.afw.geom.ellipses.Axes(
                    (match.first.getD(keyA)*lsst.afw.geom.arcseconds).asDegrees(),
                    (match.first.getD(keyB)*lsst.afw.geom.arcseconds).asDegrees(),
                    match.first.getAngle(keyTheta).asRadians()
                    )
                transform = wcs.linearizeSkyToPixel(match.first.getCoord())
                ellipse2 = lsst.afw.geom.ellipses.Quadrupole(ellipse1.transform(transform.getLinear()))
                outRecord = outCat.addNew()
                outRecord.setMomentsD(self.refEllipseKey, ellipse2)
                outRecord.setD(self.refFluxKey, calib.getFlux(match.first.getD(keyMag)))
                outRecord.setD(self.refSIndexKey, match.first.getD(keySIndex))
                prepRecord(outRecord, match.second)
        else:
            starKey = srcCat.getSchema().find("calib.psf.candidate").key
            for srcRecord in srcCat:
                if srcRecord.getFlag(starKey):
                    continue
                if where is not None and not where(src=srcRecord):
                    continue
                outRecord = outCat.addNew()
                prepRecord(outRecord, srcRecord)
        return outCat

    def fillCatalog(self, exposure, outCat):
        """For each empty ModelFitRecord in catalog, call processObject()
        """
        for n, record in enumerate(outCat):
            if n % 100 == 0:
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
