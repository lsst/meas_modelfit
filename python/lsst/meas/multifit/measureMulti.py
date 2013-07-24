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

from lsst.pipe.base import CmdLineTask, Struct, TaskError
import lsst.pex.config as pexConfig

from lsst.meas.extensions.multiShapelet import FitPsfAlgorithm
from lsst.meas.multifit import EpochImage

__all__ = ("MeasureMultiConfig", "MeasureMultiTask")

class MeasureMultiConfig(lsst.pex.config.Config):
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
        default="flat",
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


class CalexpInfo(object):
    """Basic information about a calexp that was used to generate a coadd
    
    This is intended to be compact enough that we can hold a complete list of these in memory at once
    (one per calexp/epoch). Thus it holds no image pixels. Once per object we iterate over these
    and use the information to determine which subregion of each calexp to unpersist.
    """
    def __init__(self, coaddInput, butler):
        """Construct a CalexpInfo from a CoaddInput
        
        @param[in] coaddInput: CoaddInput for a calexp
        @param[in] mapper: butler mapper containing the function getDataId(visit, ccdId)
        
        Attributes include:
        - bbox: parent bounding box (??of the whole exposure, or just the portion that overlaps the coadd??)
        - dataId: data ID dictionary
        - goodPix: ??the number of good pixels that were used in the coadd??
        - wcs: WCS
        
        coaddInput is assumed to contain at least these entries:
           (Field['PointI'](name="bbox.min", doc="bbox minimum point", units="pixels"), Key<PointI>(offset=8, nElements=2)),
           (Field['PointI'](name="bbox.max", doc="bbox maximum point", units="pixels"), Key<PointI>(offset=16, nElements=2)),
           (Field['I'](name="ccd", doc="cameraGeom CCD serial number"), Key<I>(offset=24, nElements=1)),
           (Field['L'](name="visit", doc="Foreign key for the visits (coaddTempExp) catalog"), Key<L>(offset=32, nElements=1)),
           (Field['I'](name="goodpix", doc="Number of good pixels in this CCD"), Key<I>(offset=40, nElements=1)),
        """
        self.bbox = afwGeom.Box2I(coaddInput["bbox.min"], coaddInput["bbox.max"])
        self.dataId = butler.mapper.getDataId(visit=coaddInput["visit"], ccdId=coaddInput["id"])
        self.goodPix = coaddInput["goodpix"]
        tinyBBox = afwGeom.Box2I(afwGeom.Point2I(0,0), afwGeom.Extent2I(1,1))
        tinyCalexp = butler.get("calexp_sub", bbox=tinyBBox, imageOrigin="LOCAL",
                                immediate=True, **self.dataId)
        self.wcs = tinyCalexp.getWcs()


class MeasureMultiTask(CmdLineTask):
    """Variant of MeasureImageTask for running multifit on the calexp that make up a coadd.
    
    The tasks are so different in implementation that no code is shared (yet).
    """
    ConfigClass = MeasureMultiConfig
    _DefaultName = "measureMulti"
    
    def __init__(self, **kwds):
        CmdLineTask.__init__(self, **kwds)
        self.dataPrefix = self.config.coaddName + "Coadd_"
        self.makeSubtask("sampler")
        self.schema = None # set later when we have a coadd catalog to copy
        self.basis = self.config.model.apply()
        self.prior = self.config.prior.apply()
        self.keys = {}

    def setSchema(self, coaddCat):
        """Construct self.schema; call once as soon as you have your first coadd catalog
        
        @param[in] coaddCat  ModelFitCatalog from coadd
        @raise RuntimeError if self.schema is not None
        """
        if self.schema is not None:
            raise RuntimeError("self.schema already set")

        self.log.info("Setting the schema")
        self.schemaMapper = lsst.afw.table.SchemaMapper(lsst.afw.table.SourceTable.makeMinimalSchema())
        self.schemaMapper.addMinimalSchema(lsst.meas.multifit.ModelFitTable.makeMinimalSchema())
        self.schema = self.schemaMapper.getOutputSchema()
    
        def mapKey(name, coaddCat):
            """Map one key from coaddCat to self.schema
            
            @param[in] name         name of key to map
            @param[in] coaddCat     SourceCatalog containing MeasureCoaddTask measurements
            """
            inputKey = coaddCat.getSchema()[keyName]
            self.keys[keyName] = self.schemaMapper.addMapping(inputKey)
        
        for name in ("source.ellipse", "source.center", "snr"):
            mapKey(name, coaddCat)
        if "ref" in coaddCat.getSchema():
            for name in ("ref.ellipse", "ref.center", "ref.sindex"):
                mapKey(name, coaddCat)

        def addKeys(prefix, doc):
            """Add <prefix>.ellipse and <prefix>.center keys to self.schema

            @param[in] prefix       key name prefix
            @param[in] doc          documentation prefix
            """
            self.keys["%s.ellipse" % prefix] = self.schema.addField("%s.ellipse" % prefix, type="MomentsD",
                                                                    doc=("%s ellipse" % doc))
            self.keys["%s.center" % prefix] = self.schema.addField("%s.center" % prefix, type="PointD",
                                                                   doc=("%s center position" % doc))
        addKeys("mean", "Posterior mean")
        addKeys("median", "Posterior median")

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
            outRecord = outCat.addNew()
            newRecord.assign(coaddRecord, self.schemaMapper)

        return outCat

    @lsst.pipe.base.timeMethod
    def processObject(self, outRecord, coadd, objective):
        """Process a single object.

        @param[in,out] outRecord    multi-fit ModelFitRecord
        @param[in] coadd            Coadd exposure
        @param[in] objective        multi-epoch objective (a MultiEpochObjective);
                                    optimizer objective functions that compute likelihood at a point

        @return a Struct containing various intermediate objects:
          - sampler: the Sampler object used to draw samples
          - psf: a shapelet.MultiShapeletFunction representation of the PSF
        """
        sourceCoaddCenter = outRecord.getPointD(self.keys["source.center"])
        psfModel = lsst.meas.extensions.multiShapelet.FitPsfAlgorithm.apply(
            self.config.psf.makeControl(), coadd.getPsf(), sourceCoaddCenter)
        psf = psfModel.asMultiShapelet()
        
        sampler = self.sampler.setup(
            exposure=coadd,
            center=sourceCoaddCenter,
            ellipse=outRecord.getMomentsD(self.keys["source.ellipse"]),
        )

        samples = sampler.run(objective)
        samples.applyPrior(self.prior)
        outRecord.setSamples(samples)
        mean = samples.interpret(samples.computeMean(), outRecord.getPointD(self.keys["source.center"]))
        outRecord.set(self.keys["mean.ellipse"], lsst.afw.geom.ellipses.Quadrupole(mean.getCore()))
        outRecord.set(self.keys["mean.center"], mean.getCenter())
        median = samples.interpret(samples.computeQuantiles(numpy.array([0.5])),
                                   outRecord.getPointD(self.keys["source.center"]))
        outRecord.set(self.keys["median.ellipse"], lsst.afw.geom.ellipses.Quadrupole(median.getCore()))
        outRecord.set(self.keys["median.center"], median.getCenter())

        return lsst.pipe.base.Struct(
            sampler=sampler,
            psf=psf,
        )

    def run(self, dataRef):
        """Process the catalog associated with the given coadd
        
        For each source in the catalog find those coadd input images that fully overlap the footprint,
        and fit the source shape. Combine the fits from all input images to form a multifit catalog.
        
        @param[in] dataRef: coadd data reference
        """
        butler = dataRef.getButler()
        
        inputs = self.readInputs(dataRef)
        coaddExposureInfo = inputs.coadd.getExposureInfo()
        if not coaddExposureInfo.hasCoaddInputs():
            raise TaskError("coadd does not have coaddInputs data")
        calexpInfoList = [CalexpInfo(coaddInput=coaddInput, butler=butler) \
                            for coaddInput in coaddExposureInfo.getCoaddInputs()]
        coaddWcs = coaddExposureInfo.getWcs()
        
        psfControl = self.config.psf.makeControl()

        outCat = self.prepCatalog(inputs.exposure, coaddCat=inputs.coaddCat)
            
        for outRecord in outCat:
            coaddFootprint = outRecord.getFootprint()
            sourceCoaddPos = outRecord.getPointD(self.keys["source.center"])
            
            sourceSkyPos = coaddWcs.pixelToSky(sourceCoaddPos)
            
            # process each calexp that partially overlaps this footprint
            epochImageList = []
            for calexpInfo in calexpInfoList:
                calexpFootprint = coaddFootprint.transform(coaddWcs, calexpInfo.wcs, calexpInfo.bbox)
                if calexpFootprint.getNpix() < 1:
                    # no overlapping pixels, so skip this calexp
                    continue
                calexpFootprintBBox = calexpFootprint.getBBox()
                
                calexp = butler.get(
                    "calexp_sub",
                    bbox=calexpFootprintBBox,
                    origin="PARENT",
                    immediate=True,
                    **calexpInfo.dataId)
                
                sourceCalexpPos = calexp.getWcs().skyToPixel(sourceSkyPos)
                if not calexpFootprintBBox.includes(sourceCalexpPos):
                    self.log.warn("calexp %s footprint bbox %s does not contain source at calexpPos=%s, coaddPos=%s; skipping" \
                        (calexpInfo.dataId, calexpFootprintBBox, sourceCalexpPos, sourceCoaddPos))
                    continue

                psfModel = FitPsfAlgorithm(psfControl, calexp.getPsf(), sourceCalexpPos)
                psf = psfModel.asMultiShapelet()

                epochImage = EpochImage(calexpFootprint, calexp, psf)
                epochImageList.append(epochImage)
            
            objective = MultiEpochObjective(
                self.config.objective.makeControl(),
                self.basis,
                coadd.getWcs(),
                sourceSkyPos,
                epochImageList,
            )
                
            self.processObject(outRecord=outRecord, coadd=inputs.coadd, objective=objective)
                
        
        self.writeOutputs(dataRef, inputs.coaddCat)
        return outCat

    def readInputs(self, dataRef):
        """Return inputs
        
        @param[in] dataRef: data reference for coadd
        
        @return an lsst.pipe.base.Struct containing:
          - coadd: coadd patch (lsst.afw.image.ExposureF); the images that make up this coadd are fit
          - coaddCat: catalog with model fits based on the coadd (lsst.meas.multifit.ModelFitCatalog)
        """
        return lsst.pipe.base.Struct(
            coadd = dataRef.get(self.dataPrefix[0:-1], immediate=True),
            coaddCat = dataRef.get(self.dataPrefix + "modelfits", immediate=True),
        )

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser
        """
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "deepCoadd",
            help="coadd data ID, e.g. --id tract=1 patch=2,2 filter=g")
        return parser
