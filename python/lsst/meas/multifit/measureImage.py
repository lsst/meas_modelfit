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
from lsst.meas.extensions.multiShapelet import FitPsfAlgorithm

from . import multifitLib
from .baseMeasure import BaseMeasureConfig, BaseMeasureTask
from .fitRegion import setupFitRegion

__all__ = ("MeasureImageConfig", "MeasureImageTask")

class MeasureImageConfig(BaseMeasureConfig):
    likelihood = lsst.pex.config.ConfigField(
        dtype=multifitLib.ProjectedLikelihood.ConfigClass,
        doc="Config for likelihood object that computes model probability at given parameters"
    )


class MeasureImageTask(BaseMeasureTask):
    """Driver class for S13-specific galaxy modeling work

    Like ProcessImageTask, MeasureImageTask is intended to be used as a base
    class with CCD and coadd derived-class specializations.  It is run
    after process[Ccd|Coadd|Eimage].py, and generates a single output
    catalog with the mapper name 'modelfits'.
    """

    ConfigClass = MeasureImageConfig
    dataPrefix = ""

    def __init__(self, **kwds):
        BaseMeasureTask.__init__(self, **kwds)

    def getPreviousConfig(self, butler):
        return butler.get(self._getConfigName(), tag=self.config.previous, immediate=True)

    def readInputs(self, dataRef):
        """Return a lsst.pipe.base.Struct containing the Exposure to fit and either a previous modelfits
        catalog (if config.doWarmStart) or the reference and source catalogs.
        """
        exposure = dataRef.get(self.dataPrefix + "calexp", immediate=True)
        dataset = self.dataPrefix + "modelfits"
        if self.config.previous is not None:
            prevCat = dataRef.get(self.dataPrefix + "modelfits", tag=self.config.previous, immediate=True)
            prevCat.setInterpreter(self.previous.fitter.interpreter)
            return lsst.pipe.base.Struct(
                prevCat=prevCat,
                exposure=exposure
                )
        else:
            return lsst.pipe.base.Struct(
                srcCat=dataRef.get(self.dataPrefix + "src", immediate=True),
                refCat=dataRef.get("refcat", immediate=True),
                exposure=exposure
                )

    def prepCatalog(self, inputs):
        """Prepare and return the output catalog, doing everything but the actual fitting.

        If config.doWarmStart, this just returns the previous modelfits catalog we're using
        for the warm start.  If not, it crease a new modelfits catalog by matching the source
        and reference catalogs, transforming their fields to the fit coordinate system and
        amplitude units.
        """
        outCat = multifitLib.ModelFitCatalog(self.makeTable())
        refCat = inputs.refCat
        srcCat = inputs.srcCat

        exposureWcs = inputs.exposure.getWcs()

        # SchemaMapper will transfer ID, Coord, Footprint (we'll overwrite the Coord with that from RefCat)
        mapper = lsst.afw.table.SchemaMapper(lsst.afw.table.SourceTable.makeMinimalSchema())
        mapper.addMinimalSchema(lsst.meas.multifit.ModelFitTable.makeMinimalSchema())

        # extract some keys from ref catalog for later use
        keyA = refCat.getSchema().find("ellipse.a").key
        keyB = refCat.getSchema().find("ellipse.b").key
        keyTheta = refCat.getSchema().find("ellipse.theta").key
        keyMag = refCat.getSchema().find("mag.%s" % inputs.exposure.getFilter().getName()).key
        keySIndex = refCat.getSchema().find("sindex").key

        # Do a spatial match between srcCat and refCat to determine what to fit: we use
        # the refCat (which has only galaxies) to avoid fitting stars.
        matches = lsst.afw.table.matchRaDec(refCat, srcCat, 1.0*lsst.afw.geom.arcseconds)
        for match in matches:
            srcRecord = match.second
            refRecord = match.first
            outRecord = outCat.addNew()

            # Start by setting some miscellaneous calculated fields
            outRecord.assign(srcRecord, mapper)
            outRecord.setD(self.keys["snr"], srcRecord.getApFlux() / srcRecord.getApFluxErr())
            outRecord.setCoord(refRecord.getCoord())
            outRecord.setPointD(self.keys["center"], exposureWcs.skyToPixel(refRecord.getCoord()))

            # Next we determine the pixel region we want to fit.
            outRecord.setFootprint(setupFitRegion(self.config.fitRegion, inputs.exposure, srcRecord))

            # Setup the coordinate and photometric systems to use for the parameters
            units = self.makeUnitSystem(outRecord, outRecord.getCoord(), refRecord.get(keyMag))

            # Now we'll transform the refCat ellipse to the parameter coordinate system.
            ellipse1 = lsst.afw.geom.ellipses.Ellipse(
                lsst.afw.geom.ellipses.Axes(
                    (refRecord.getD(keyA)*lsst.afw.geom.arcseconds).asDegrees(),
                    (refRecord.getD(keyB)*lsst.afw.geom.arcseconds).asDegrees(),
                    (refRecord.getD(keyTheta)*lsst.afw.geom.degrees).asRadians() - 0.5*numpy.pi
                    ),
                refRecord.getCoord().getPosition(lsst.afw.geom.degrees)
                )
            transform = units.wcs.linearizeSkyToPixel(refRecord.getCoord())
            ellipse2 = ellipse1.transform(transform)

            # We now transform this ellipse and the refCat fluxes into the parameters defined by the model.
            ellipses = lsst.meas.multifit.Model.EllipseVector()
            ellipses.append(ellipse2)
            nonlinear = outRecord[self.keys["initial.nonlinear"]]
            self.model.readEllipses(ellipses, nonlinear, outRecord[self.keys["fixed"]])

            # this flux->amplitudes conversion assumes the ref catalog is single-component, and that the
            # first component of the model is what that corresponds to; we may need to generalize this
            amplitudes = outRecord[self.keys["initial.amplitudes"]]
            amplitudes[:] = 0.0
            amplitudes[0] = 1.0

            # Finally, we tell the fitter to initialize the record, allowing it to do any fitter-specific
            # bootstrapping of the record.
            self.fitter.initialize(outRecord)
        return outCat

    def makeLikelihood(self, inputs, record):
        """Create a Likelihood object for a single object.

        The MeasureImage implementation creates a ProjectedLikelihood with data from a single
        exposure.
        """
        psfModel = FitPsfAlgorithm.apply(self.config.psf.makeControl(), inputs.exposure.getPsf(),
                                         record.get(self.keys["center"]))
        psf = psfModel.asMultiShapelet()
        return multifitLib.ProjectedLikelihood(
            self.model, record[self.keys["fixed"]],
            self.getUnitSystem(record),
            record.getCoord(),
            inputs.exposure,
            record.getFootprint(),
            psf,
            self.config.likelihood.makeControl()
            )

    def writeOutputs(self, dataRef, outCat):
        dataRef.put(outCat, self.dataPrefix + "modelfits", tag=self.config.tag)

    def getSchemaCatalogs(self):
        """Return a dict of empty catalogs for each catalog dataset produced by this task."""
        return {self.dataPrefix + "modelfits": multifitLib.ModelFitCatalog(self.makeTable())}
