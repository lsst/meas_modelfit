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

from .samplers import AdaptiveImportanceSamplerTask
from .fitRegion import setupFitRegion
from .models import *
from .priors import *

__all__ = ("BaseMeasureConfig", "BaseMeasureTask", "InitialMeasureConfig", "InitialMeasureTask")

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

    def setDefaults(self):
        self.psf.innerOrder = 4
        self.psf.outerOrder = 0

class BaseMeasureTask(lsst.pipe.base.CmdLineTask):
    """Base class for MeasureImageTask and MeasureMultiTask to aggregate shared code"""

    def __init__(self, **kwds):
        lsst.pipe.base.CmdLineTask.__init__(self, **kwds)
        self.keys = {}
        self.schema = None

    def addEllipseFields(self, prefix, doc):
        """Add <prefix>.ellipse and <prefix>.center keys to self.schema, saving keys in self.keys
        @param[in] prefix       key name prefix
        @param[in] doc          documentation prefix
        """
        self.keys["%s.ellipse" % prefix] = self.schema.addField("%s.ellipse" % prefix, type="MomentsD",
                                                                doc=("%s ellipse" % doc))
        self.keys["%s.center" % prefix] = self.schema.addField("%s.center" % prefix, type="PointD",
                                                               doc=("%s center position" % doc))
    def addMeasuredFields(self):
        """Add fields to self.schema for quantities derived from the SampleSet
        """
        self.addEllipseFields("mean", "Posterior mean")
        self.addEllipseFields("median", "Posterior median")

    def fillMeasuredFields(self, record):
        samples = record.getSamples()
        parameterDef = samples.getParameterDefinition()
        mean = parameterDef.makeEllipse(samples.computeMean(), record.getPointD(self.keys["source.center"]))
        record.set(self.keys["mean.ellipse"], lsst.afw.geom.ellipses.Quadrupole(mean.getCore()))
        record.set(self.keys["mean.center"], mean.getCenter())
        median = parameterDef.makeEllipse(samples.computeQuantiles(numpy.array([0.5])),
                                          record.getPointD(self.keys["source.center"]))
        record.set(self.keys["median.ellipse"], lsst.afw.geom.ellipses.Quadrupole(median.getCore()))
        record.set(self.keys["median.center"], median.getCenter())

class InitialMeasureConfig(BaseMeasureConfig):
    fitRegion = lsst.pex.config.ConfigField(
        dtype=setupFitRegion.ConfigClass,
        doc="Parameters that control which pixels to include in the model fit"
    )

class InitialMeasureTask(BaseMeasureTask):
    """Intermediate base class for measurement tasks that start by matching sources and references.
    """

    def __init__(self, **kwds):
        BaseMeasureTask.__init__(self, **kwds)
        self.keys = {}
        self.refKeys = {}
        self.schema = None
        self.srcMapper = None

    def setupInitialSchema(self, exposure, srcCat, refCat):
        """Attach Schema and SchemaMapper attributes to self (.schema and self.srcMapper, respectively)
        for use with prepInitialCatalog().
        """
        srcSchema = srcCat.getSchema()
        refSchema = refCat.getSchema()
        self.srcMapper = lsst.afw.table.SchemaMapper(srcSchema)
        self.srcMapper.addMinimalSchema(lsst.meas.multifit.ModelFitTable.makeMinimalSchema()) # just ID
        self.srcMapper.addMapping(srcCat.getTable().getShapeKey(), "source.ellipse")
        self.srcMapper.addMapping(srcCat.getTable().getCentroidKey(), "source.center")
        self.schema = self.srcMapper.getOutputSchema()
        self.keys["snr"] = self.schema.addField("snr", type=float,
                                                doc="signal to noise ratio from source apFlux/apFluxErr")
        self.addEllipseFields("ref", "Reference catalog")
        self.keys["ref.sindex"] = self.schema.addField("ref.sindex", type=float,
                                                       doc="Reference catalog Sersic index")
        self.keys["ref.flux"] = self.schema.addField("ref.flux", type=float,
                                                     doc="Reference catalog flux")
        self.refKeys["ellipse.a"] = refSchema.find("ellipse.a").key
        self.refKeys["ellipse.b"] = refSchema.find("ellipse.b").key
        self.refKeys["ellipse.theta"] = refSchema.find("ellipse.theta").key
        self.refKeys["mag"] = refSchema.find("mag.%s" % exposure.getFilter().getName()).key
        self.refKeys["sindex"] = refSchema.find("sindex").key
        self.addMeasuredFields()

    def prepInitialCatalog(self, exposure, srcCat, refCat, where=None):
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
        if self.schema is None:
            self.setupInitialSchema(exposure=exposure, srcSchema=srcCat.getSchema(),
                                    refSchema=refCat.getSchema())
        outCat = ModelFitCatalog(self.schema)
        matches = lsst.afw.table.matchRaDec(refCat, srcCat, 1.0*lsst.afw.geom.arcseconds)
        wcs = exposure.getWcs()
        calib = exposure.getCalib()
        for match in matches:
            refRecord = match.first
            srcRecord = match.second
            if where is not None and not where(src=srcRecord, ref=refRecord):
                continue
            if srcRecord.getShapeFlag():
                continue
            ellipse1 = lsst.afw.geom.ellipses.Axes(
                (refRecord.getD(self.refKeys["ellipse.a"])*lsst.afw.geom.arcseconds).asDegrees(),
                (refRecord.getD(self.refKeys["ellipse.b"])*lsst.afw.geom.arcseconds).asDegrees(),
                (refRecord.getD(self.refKeys["ellipse.theta"])*lsst.afw.geom.degrees).asRadians()
                - 0.5*numpy.pi
                )
            transform = wcs.linearizeSkyToPixel(refRecord.getCoord())
            ellipse2 = lsst.afw.geom.ellipses.Quadrupole(ellipse1.transform(transform.getLinear()))
            outRecord = outCat.addNew()
            outRecord.assign(srcRecord, self.srcMapper)
            outRecord.setMomentsD(self.keys["ref.ellipse"], ellipse2)
            outRecord.setD(self.keys["ref.flux"], calib.getFlux(refRecord.getD(self.refKeys["mag"])))
            outRecord.setD(self.keys["ref.sindex"], refRecord.getD(self.refKeys["sindex"]))
            outRecord.setPointD(self.keys["ref.center"], wcs.skyToPixel(refRecord.getCoord()))
            outRecord.setFootprint(setupFitRegion(self.config.fitRegion, exposure, srcRecord))
            outRecord.setD(self.keys["snr"], srcRecord.getApFlux() / srcRecord.getApFluxErr())
        return outCat
