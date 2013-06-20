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

from .measureImage import MeasureImageTask

__all__ = ("MeasureMultiConfig", "MeasureMultiTask")

class MeasureMultiConfig(pexConfig.Config):
    coaddName = pexConfig.Field(
        doc = "Coadd name: typically one of deep or goodSeeing.",
        dtype = str,
        default = "deep",
    )
    sampler = lsst.pex.config.ConfigurableField(
        target=NaiveGridSamplerTask,
        doc="Subtask that generates samples from the probability of a galaxy model given image data"
    )
    objective = lsst.pex.config.ConfigField(
        dtype=SingleEpochObjective.ConfigClass,
        doc="Config for objective object that computes model probability at given parameters"
    )
    useRefCat = lsst.pex.config.Field(
        dtype=bool,
        default=True,
        doc="Whether to use the reference catalog to identify objects to fit"
    )


class CalexpInfo(object):
    """Basic information about a calexp that was used to generate a coadd
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
    """Specialization of MeasureImageTask for running multifit on the calexp that make up a coadd.
    """
    ConfigClass = MeasureMultiConfig
    _DefaultName = "measureMulti"
    
    def __init__(self, **kwds):
        CmdLineTask.__init__(self, **kwds)
        self.dataPrefix = self.config.coaddName + "Coadd_"

    @lsst.pipe.base.timeMethod
    def processObject(self, exposureInfo, record):
        """Process a single object.

        !!!COPIED FROM MeasureImageTask; NEEDS MAJOR CHANGES!!!
        
        @param[in] exposureInfo ExposureInfo for the coadd
        @param[in,out] record   ModelFitRecord to fill

        @return a Struct containing:
        - objective: the Objective object used to do the fitting
        - ?
        """
        if not exposure.hasCoaddInputs():
            raise pipeBase.TaskError("exposure does not have coaddInputs data")
        
        coaddInputs = exposure.getCoaddInputs()
        ccdCat = coaddInputs.ccds
            
        psfModel = lsst.meas.extensions.multiShapelet.FitPsfModel(self.config.psf.makeControl(), record)
        psf = psfModel.asMultiShapelet()
        sampler = self.sampler.setup(exposure=exposure, center=record.getPointD(self.keys["source.center"]),
                                     ellipse=record.getMomentsD(self.keys["source.ellipse"]))
        objective = SingleEpochObjective(
            self.config.objective.makeControl(), self.basis, psf,
            exposure.getMaskedImage(), record.getFootprint()
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

        return lsst.pipe.base.Struct(
            objective=objective,
            sampler=sampler,
            psf=psf,
            record=record,
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
            
        for i, modelFit in enumerate(inputs.modelFitCat):
            objCoaddBBoxD = afwGeom.Box2D(modelFit.getFootprint().getBBox())
            objCoordList = [coaddWcs.pixelToSky(pixPos) for pixPos in objCoaddBBoxD.getCorners()]
            
            # process each calexp that fully overlaps this footprint
            calexpList = []
            for calexpInfo in calexpInfoList:
                objCalexpBBox = [calexpInfo.wcs.skyToPixel(skyPos) for skyPos in objCoordList]
                if not calexpInfo.bbox.includes(objCalexpBBox):
                    # object is not fully on this calexp, so skip the calexp
                    continue

                # ???do we just need the objCalexpBBox pixels, or should we grow it to provide some margin;
                # and if we grow it, do we grow before testing for inclusion or afterwards???
                calexp = butler.get("calexp_sub", bbox=objCalexpBBox, origin="PARENT", **calexpInfo.dataId)
                calexpList.append(calexp)
                
            self.processObject(modelFit, calexpList, ?????)
                
        
        ???self.writeOutputs(dataRef, outCat)???
        return outCat

    def readInputs(self, dataRef):
        """Return inputs
        
        @param[in] dataRef: data reference for coadd
        
        @return an lsst.pipe.base.Struct containing:
          - coadd: coadd patch (lsst.afw.image.ExposureF); the images that make up this coadd are fit
          - modelFitCat: catalog with model fits based on the coadd (lsst.meas.multifit.ModelFitCatalog)
          - refCat: catalog with with truth values (lsst.afw.table.SimpleCatalog);
                may be or None if not config.useRefCat
        """
        try:
            refCat = dataRef.get("refcat", immediate=True)
        except:
            refCat = None
        return lsst.pipe.base.Struct(
            coadd = dataRef.get(self.dataPrefix[0:-1], immediate=True),
            modelFitCat = dataRef.get(self.dataPrefix + "modelfits", immediate=True),
            refCat = refCat,
        )

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser
        """
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "deepCoadd",
            help="coadd data ID, e.g. --id tract=1 patch=2,2 filter=g")
        return parser
