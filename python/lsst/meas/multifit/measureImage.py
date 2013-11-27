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
    doUseDataUnitSystem = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        doc=("If True, use the WCS and Calib of the Exposure being fit to define the units for the "
             "parameters to be fit.  The default is False, which sets up a local unit system that "
             "ensures all parameters are of order unity and the prior unit system matches exactly.")
    )
    doWarmStart = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        doc=("If True, load a previous modelfits catalog and use its attached proposal distributions "
             "as the initial proposal distributions, instead of starting from a match of the source "
             "and reference catalogs.  NOTE: Also causes fit region footprints to be based on the previous "
             "modelfits footprints, instead of the original detection footprints")
    )

    def makeUnitSystem(self, position, magnitude, exposure):
        """Create the UnitSystem that defines the geometric and photometric units for the parameters.

        All arguments are always passed, but which are used depends on the value of doUseDataUnitSystem.

        @param[in] position    Position of the object to be fit, as an afw::coord::Coord, used to create
                               a local TAN WCS at the position of the object with arcsecond pixels.  Used
                               only if doUseDataUnitSystem is False.
        @param[in] magnitude   Approximate magnitude of the object to be fit; used to scale the flux units
                               such that amplitude parameters are of order 1.  Used only if
                               doUseDataUnitSystem is False.
        @param[in] exposure    The exposure being providing data for the fit.  Used only if
                               doUseDataUnitSystem is True.

        """
        if self.doUseDataUnitSystem:
            return multifitLib.UnitSystem(exposure)
        else:
            return multifitLib.UnitSystem(position, magnitude)

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
        if not self.config.doUseDataUnitSystem:
            self.keys["sys.position"] = self.schema.addField(
                "sys.position", type="Coord",
                doc="nominal position used to construct UnitSystem for parameters"
                )
            self.keys["sys.magnitude"] = self.schema.addField(
                "sys.mag", type=float,
                doc="nominal magnitude used to construct UnitSystem for parameters"
                )

    def readInputs(self, dataRef):
        """Return a lsst.pipe.base.Struct containing the Exposure to fit and either a previous modelfits
        catalog (if config.doWarmStart) or the reference and source catalogs.
        """
        exposure = dataRef.get(self.dataPrefix + "calexp", immediate=True)
        if self.config.doWarmStart:
            # TODO: load previous config, makes sure it's compatible (in particular, has same unit system)
            return lsst.pipe.base.Struct(
                prevCat=dataRef.get(self.dataPrefix + "modelfits", immediate=True),
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

        prior = self.config.makePrior(exposure=inputs.exposure)
        interpreter = self.fitter.interpreter.clone()
        interpreter.setPrior(prior)

        if self.config.doWarmStart:
            # Interpreters aren't persisted with the ModelFitCatalog, since they can be
            # reconstructed entirely from Config, so we reattach one here.
            inputs.prevCat.table.setInterpreter(interpreter)
            return inputs.prevCat

        table = self.makeTable()
        table.setInterpreter(interpreter)
        outCat = multifitLib.ModelFitCatalog(table)
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
            outRecord.setPointD(self.keys["ref.center"], exposureWcs.skyToPixel(refRecord.getCoord()))

            # Next we determine the pixel region we want to fit.
            outRecord.setFootprint(setupFitRegion(self.config.fitRegion, inputs.exposure, srcRecord))

            # Setup the coordinate and photometric systems to use for the fit
            fitSys = self.config.makeFitSystem(outRecord.getCoord(), exposure=inputs.exposure)

            # Now we'll transform the refCat ellipse to the fit coordinate system.
            ellipse1 = lsst.afw.geom.ellipses.Ellipse(
                lsst.afw.geom.ellipses.Axes(
                    (refRecord.getD(keyA)*lsst.afw.geom.arcseconds).asDegrees(),
                    (refRecord.getD(keyB)*lsst.afw.geom.arcseconds).asDegrees(),
                    (refRecord.getD(keyTheta)*lsst.afw.geom.degrees).asRadians() - 0.5*numpy.pi
                    ),
                refRecord.getCoord().getPosition(lsst.afw.geom.degrees)
                )
            transform = fitSys.wcs.linearizeSkyToPixel(refRecord.getCoord())
            ellipse2 = ellipse1.transform(transform)

            # We now transform this ellipse and the refCat fluxes into the parameters defined by the model.
            ellipses = lsst.meas.multifit.Model.EllipseVector()
            ellipses.append(ellipse2)
            self.model.readEllipses(
                ellipses,
                outRecord[self.keys["ref.nonlinear"]],
                outRecord[self.keys["ref.fixed"]]
                )

            # this flux->amplitudes conversion assumes the ref catalog is single-component, and that the
            # first component of the model is what that corresponds to; we may need to generalize this
            flux = fitSys.calib.getFlux(refRecord.get(keyMag))
            amplitudes = outRecord[self.keys["ref.amplitudes"]]
            amplitudes[:] = 0.0
            amplitudes[0] = flux

            # Finally, we tell the fitter to initialize the record, which includes setting "ref.parameters",
            # as only the fitter knows what parameters it will actually fit and how those map to
            # "ref.nonlinear" and "ref.amplitudes".  Depending on the fitter, it will probably attach
            # an initial PDF too.
            self.fitter.initialize(outRecord)
        return outCat

    def makeLikelihood(self, inputs, record):
        """Create a Likelihood object for a single object.

        The MeasureImage implementation creates a ProjectedLikelihood with data from a single
        exposure.
        """
        psfModel = FitPsfAlgorithm.apply(self.config.psf.makeControl(), inputs.exposure.getPsf(),
                                         record.get(self.keys["ref.center"]))
        psf = psfModel.asMultiShapelet()
        return multifitLib.ProjectedLikelihood(
            self.model, record[self.keys["ref.fixed"]],
            self.config.makeFitSystem(record.getCoord(), exposure=inputs.exposure),
            record.getCoord(),
            inputs.exposure,
            record.getFootprint(),
            psf,
            self.config.likelihood.makeControl()
            )

    def writeOutputs(self, dataRef, outCat):
        dataRef.put(outCat, self.dataPrefix + "modelfits")

    def getSchemaCatalogs(self):
        """Return a dict of empty catalogs for each catalog dataset produced by this task."""
        return {self.dataPrefix + "modelfits": multifitLib.ModelFitCatalog(self.makeTable())}
