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

from . import multifitLib
from .baseMeasure import BaseMeasureConfig, BaseMeasureTask

try:
    from lsst.meas.mosaic import applyMosaicResults
except ImportError:
    applyMosaicResults = None

__all__ = ("MeasureMultiConfig", "MeasureMultiTask")

class MeasureMultiConfig(BaseMeasureConfig):
    coaddName = lsst.pex.config.Field(
        doc="coadd name: typically one of deep or goodSeeing",
        dtype=str,
        default="deep",
    )
    likelihood = lsst.pex.config.ConfigField(
        dtype=multifitLib.ProjectedLikelihood.ConfigClass,
        doc="Config for likelihood object that computes model probability at given parameters"
    )
    minPixels = lsst.pex.config.Field(
        doc="minimum number of pixels in a calexp footprint to use that calexp for a given galaxy",
        dtype=int,
        default=5,
    )
    usePreviousMultiFit = lsst.pex.config.Field(
        dtype=bool,
        doc="If True, do a warm-start from a previous MeasureMulti run; if False, use MeasureCoadd outputs",
        default=False
    )
    doApplyUberCal = lsst.pex.config.Field(
        dtype=bool,
        doc="Apply meas_mosaic ubercal results to input calexps?",
        default=True
    )

class MeasureMultiTask(BaseMeasureTask):
    """Variant of BaseMeasureTask for running multifit on the calexps that make up a coadd.

    The tasks are so different in implementation that no code is shared (yet).
    """
    ConfigClass = MeasureMultiConfig
    _DefaultName = "measureMulti"

    def __init__(self, **kwds):
        BaseMeasureTask.__init__(self, **kwds)
        self.outputName = self.config.coaddName + "Multi_modelfits"

    def readInputs(self, dataRef):
        """Return inputs, and attach to the task a closure method (readInputExposure)
        that loads a calexp subimage given an ExposureRecord and a bounding box.

        @param[in] dataRef: data reference for coadd

        @return an lsst.pipe.base.Struct containing:
          - prevCat: ModelFitCatalog used for a "warm start" for the fitting
          - exposureCat: catalog of ExposureRecords that determines which calexps can be used in the fitting
          - footprintWcs: Wcs of the Footprints attached to prevCat records
          - readInputExposure: a closure method, used to load individual calexp subimages
        """
        if self.config.usePreviousMultiFit:
            prevCat = dataRef.get(self.outputName, immediate=True)
        else:
            prevCat = dataRef.get(self.config.coaddName + "Coadd_modelfits", immediate=True)
        # TODO: when possible, just load the coaddInputCat and Wcs, not the full coadd
        coadd = dataRef.get(self.config.coaddName, immediate=True)
        footprintWcs = coadd.getWcs()
        exposureCat = coadd.getInfo().getCoaddInputs().ccds
        exposureSchema = coaddInputCat.getSchema()
        visitKey = exposureSchema.find("visit").key
        ccdKey = exposureSchema.find("ccd").key
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
        return lsst.pipe.base.Struct(
            prevCat=prevCat,
            exposureCat=exposureCat,
            footprintWcs=footprintWcs,
            readInputExposure=readInputExposure
        )

    def prepCatalog(self, inputs):
        """Prepare and return the output catalog, doing everything but the actual fitting.

        After this step, each output record should be in a state such that makeLikelihood and
        fitter.run() may be called on it.
        """
        return inputs.prevCat

    @lsst.pipe.base.timeMethod
    def makeLikelihood(self, inputs, record):

        # process each calexp that partially overlaps this footprint
        epochFootprintList = multifitLib.EpochFootprintVector()

        psfCtrl = self.config.psf.makeControl()
        fitWcs = self.makeFitWcs(record.getCoord())

        for exposureRecord in inputs.exposureCat:
            calexpFootprint = coaddFootprint.transform(inputs.footprintWcs, exposureRecord.getWcs(),
                                                       exposureRecord.getBBox())

            if calexpFootprint.getArea() < self.config.minPixels:
                continue
            calexpFootprintBBox = calexpFootprint.getBBox()
            assert not calexpFootprintBBox.isEmpty() # verify that #2979 is fixed in afw

            calexp = self.readInputExposure(record=exposureRecord, bbox=calexpFootprintBBox)

            sourceCalexpPos = calexp.getWcs().skyToPixel(record.getCoord())

            psfModel = FitPsfAlgorithm.apply(psfCtrl, calexp.getPsf(), sourceCalexpPos)
            psf = psfModel.asMultiShapelet()

            epochFootprint = multifitLib.EpochFootprint(calexpFootprint, calexp, psf)
            epochFootprintList.append(epochFootprint)

        return MultiEpochLikelihood(
            self.model, record.get(self.keys["fixed"]),
            fitWcs, self.fitCalib,
            record.getCoord(),
            epochFootprintList,
            self.config.likelihood.makeControl()
        )

    def writeOutputs(self, dataRef, outCat):
        dataRef.put(outCat, self.outputName)

    @classmethod
    def _makeArgumentParser(cls):
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

    def getSchemaCatalogs(self):
        """Return a dict of empty catalogs for each catalog dataset produced by this task."""
        return {self.outputName: ModelFitCatalog(self.makeTable())}
