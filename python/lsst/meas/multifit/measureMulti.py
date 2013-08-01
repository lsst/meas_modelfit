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

from lsst.pipe.base import CmdLineTask, Struct, TaskError
import lsst.pex.config
import lsst.afw.table # temporary hack; see PSF HACK below
import lsst.afw.geom as afwGeom
from lsst.meas.extensions.multiShapelet import FitPsfAlgorithm
from .multifitLib import VectorEpochFootprint, EpochFootprint, MultiEpochObjective, ModelFitCatalog, \
    ModelFitTable
from .measureImage import BaseMeasureConfig 

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

class CalexpInfo(object):
    """Basic information about a calexp that was used to generate a coadd
    
    This is intended to be compact enough that we can hold a complete list of these in memory at once
    (one per calexp/epoch). Thus it holds no image pixels. Once per object we iterate over these
    and use the information to determine which subregion of each calexp to unpersist.
    """
    @classmethod
    def setKeys(cls, schema):
        """Obtain the required keys from the specified schema.
        
        Must be called before instantiating any instances of this class.
        """
        cls.keys = dict()
        for keyName in ("bbox.min", "bbox.max", "visit", "ccd", "goodpix"):
            cls.keys[keyName] = schema.find(keyName).getKey()
        
    def __init__(self, coaddInput, butler):
        """Construct a CalexpInfo from a CoaddInput
        
        @param[in] coaddInput: CoaddInput for a calexp
        @param[in] butler: data butler whose mapper contains the function getDataId(visit, ccdId)
        
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
        self.bbox = afwGeom.Box2I(coaddInput.get(self.keys["bbox.min"]), coaddInput.get(self.keys["bbox.max"]))
        self.dataId = butler.mapper.getDataId(
            visit=coaddInput.get(self.keys["visit"]),
            ccdId=coaddInput.get(self.keys["ccd"]),
        )
        self.goodPix = coaddInput[self.keys["goodpix"]]
        tinyBBox = afwGeom.Box2I(afwGeom.Point2I(0,0), afwGeom.Extent2I(1,1))
        tinyCalexp = butler.get(
            "calexp_sub",
            bbox=tinyBBox,
            imageOrigin="LOCAL",
            immediate=True,
            **self.dataId
        )
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
        self.psfControl = self.config.psf.makeControl()
        self.keys = {}

    def run(self, dataRef):
        """Process the catalog associated with the given coadd
        
        For each source in the catalog find those coadd input images that fully overlap the footprint,
        and fit the source shape. Combine the fits from all input images to form a multifit catalog.
        
        @param[in] dataRef: coadd data reference
        """
        inputs = self.readInputs(dataRef)
        outCat = self.prepCatalog(coaddCat=inputs.coaddCat)

        if not self.config.prepOnly:
            butler = dataRef.getButler()
            calexpInfoList = self.getCalexpInfoList(butler=butler, coadd=inputs.coadd)
            for record in outCat:
                self.processObject(record=record, butler=butler, coadd=inputs.coadd, calexpInfoList=calexpInfoList)
                
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

        return outCat

    def getCalexpInfoList(self, butler, coadd):
        """Get a list of CalexpInfo objects, one per calexp in the coadd
        
        @param[in] butler: data butler
        @param[in] coadd: coadd exposure
        """
        coaddExposureInfo = coadd.getInfo()
        if not coaddExposureInfo.hasCoaddInputs():
            raise TaskError("coadd exposure info does not contain coadd inputs data")
        coaddInputs = coaddExposureInfo.getCoaddInputs()
        CalexpInfo.setKeys(coaddInputs.ccds.schema)
        return [CalexpInfo(coaddInput=coaddInput, butler=butler) for coaddInput in coaddInputs.ccds]

    @lsst.pipe.base.timeMethod
    def processObject(self, record, butler, coadd, calexpInfoList):
        """Process a single object.

        @param[in,out] record   multi-fit ModelFitRecord
        @param[in] butler       data butler
        @param[in] coadd        Coadd exposure
        @param[in] objective    multi-epoch objective (a MultiEpochObjective);
                                optimizer objective functions that compute likelihood at a point

        @return a Struct containing various intermediate objects:
          - objective   the Objective object used to evaluate likelihoods
          - sampler     the Sampler object used to draw samples
          - record      the output record (identical to the record argument, which is modified in-place)
        """
        objective = self.makeObjective(
            record=record,
            butler=butler,
            coadd=coadd,
            calexpInfoList=calexpInfoList,
        )

        sourceCoaddCenter = record.getPointD(self.keys["source.center"])
        
        sampler = self.sampler.setup(
            exposure=coadd,
            center=sourceCoaddCenter,
            ellipse=record.getMomentsD(self.keys["source.ellipse"]),
        )

        samples = sampler.run(objective)
        samples.applyPrior(self.prior)
        record.setSamples(samples)
        mean = samples.interpret(samples.computeMean(), record.getPointD(self.keys["source.center"]))
        record.set(self.keys["mean.ellipse"], lsst.afw.geom.ellipses.Quadrupole(mean.getCore()))
        record.set(self.keys["mean.center"], mean.getCenter())
        median = samples.interpret(samples.computeQuantiles(numpy.array([0.5])),
                                   record.getPointD(self.keys["source.center"]))
        record.set(self.keys["median.ellipse"], lsst.afw.geom.ellipses.Quadrupole(median.getCore()))
        record.set(self.keys["median.center"], median.getCenter())

        return lsst.pipe.base.Struct(objective=objective, sampler=sampler, record=record)

    @lsst.pipe.base.timeMethod
    def makeObjective(self, record, butler, coadd, calexpInfoList):
        """Construct a MultiEpochObjective from the calexp footprints

        @param[in] record   multi-fit ModelFitRecord
        @param[in] butler       data butler
        @param[in] coadd        Coadd exposure
        @param[in] calexpInfoList    list of CalexpInfo, one per calexp in the coadd
        """
        coaddFootprint = record.getFootprint()
        sourceCoaddPos = record.getPointD(self.keys["source.center"])
        
        coaddWcs = coadd.getWcs()
        sourceSkyPos = coaddWcs.pixelToSky(sourceCoaddPos)

        # process each calexp that partially overlaps this footprint
        epochFootprintList = VectorEpochFootprint()
        for calexpInfo in calexpInfoList:
            calexpFootprint = coaddFootprint.transform(coaddWcs, calexpInfo.wcs, calexpInfo.bbox)
            calexpFootprint.clipTo(calexpInfo.bbox) # Footprint.transform does not clip
            
            # due to ticket #2979 this test is invalid -- getNpix includes clipped pixels!
            # so also test that the footprint's bbox is not empty
            if calexpFootprint.getNpix() < self.config.minPixels:
                # no overlapping pixels, so skip this calexp
                continue
            calexpFootprintBBox = calexpFootprint.getBBox()
            if calexpFootprintBBox.isEmpty(): # temporary hack due to ticket #2979
                continue
            
            calexp = butler.get(
                "calexp_sub",
                bbox=calexpFootprintBBox,
                origin="PARENT",
                immediate=True,
                **calexpInfo.dataId
            )
            
            sourceCalexpPos = calexp.getWcs().skyToPixel(sourceSkyPos)
# the following test is apparently not necessary; we just want a few valid pixels            
#             if not calexpFootprintBBox.contains(afwGeom.Point2I(sourceCalexpPos)):
#                 self.log.warn("calexp %s footprint bbox %s does not contain source at calexpPos=%s, coaddPos=%s; skipping" % \
#                     (calexpInfo.dataId, calexpFootprintBBox, sourceCalexpPos, sourceCoaddPos))
#                 continue

            # PSF HACK the next line should work; since it doesn't, use the mess below it instead
            #psfModel = FitPsfAlgorithm.apply(self.psfControl, calexp.getPsf(), sourceCalexpPos)
            schema = lsst.afw.table.Schema()
            psfFitter = FitPsfAlgorithm(self.psfControl, schema)
            tempTable = lsst.afw.table.BaseTable.make(schema)
            tempRecord = tempTable.makeRecord()
            psfModel = psfFitter.apply(tempRecord, calexp.getPsf(), sourceCalexpPos)

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
        
        for name in ("source.ellipse", "source.center", "snr"):
            mapKey(name)
        if "ref" in coaddCat.getSchema():
            for name in ("ref.ellipse", "ref.center", "ref.sindex"):
                mapKey(name)

        self.schema = self.schemaMapper.getOutputSchema()

        def addKeys(prefix, doc):
            """Add <prefix>.ellipse and <prefix>.center keys to self.schemaMapper
            @param[in] prefix       key name prefix
            @param[in] doc          documentation prefix
            """
            self.keys["%s.ellipse" % prefix] = self.schema.addField("%s.ellipse" % prefix, type="MomentsD",
                                                                    doc=("%s ellipse" % doc))
            self.keys["%s.center" % prefix] = self.schema.addField("%s.center" % prefix, type="PointD",
                                                                   doc=("%s center position" % doc))
        addKeys("mean", "Posterior mean")
        addKeys("median", "Posterior median")

    def writeOutputs(self, dataRef, outCat):
        """Write task outputs using the butler.
        """
        self.log.info("Writing output catalog")
        outCat.sort()  # want to sort by ID before saving, so when we load it's contiguous
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
