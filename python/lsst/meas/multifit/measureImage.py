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

from . import multifitLib
from .baseMeasure import BaseMeasureConfig, BaseMeasureTask

__all__ = ("MeasureImageConfig", "MeasureImageTask")

class MeasureImageConfig(BaseMeasureConfig):
    likelihood = lsst.pex.config.ConfigField(
        dtype=multifitLib.ProjectedLikelihood.ConfigClass,
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

class MeasureImageTask(BaseMeasureTask):
    """Driver class for S13-specific galaxy modeling work

    Like ProcessImageTask, MeasureImageTask is intended to be used as a base
    class with CCD and coadd derived-class specializations.  It is run
    after process[Ccd|Coadd|Eimage].py, and generates a single output
    catalog with the mapper name 'modelfits'.
    """

    ConfigClass = MeasureImageConfig
    dataPrefix = ""

    def readInputs(self, dataRef):
        """Return a lsst.pipe.base.Struct containing the Exposure to fit and either a previous modelfits
        catalog (if config.doWarmStart) or the reference and source catalogs.
        """
        exposure = dataRef.get(self.dataPrefix + "calexp", immediate=True),
        if self.config.doWarmStart:
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
        if self.config.doWarmStart:
            return inputs.prevCat

        outCat = multifitLib.ModelFitCatalog(self.makeTable())
        refCat = inputs.refCat
        srcCat = inputs.srcCat

        # SchemaMapper will transfer ID, Coord, Footprint (we'll overwrite the Coord with that from RefCat)
        mapper = lsst.afw.table.SchemaMapper(lsst.afw.table.SourceTable.makeMinimalSchema())
        mapper.addMinimalSchema(lsst.meas.multifit.ModelFitTable.makeMinimalSchema())

        # extract some keys from ref catalog for later use
        keyA = refCat.getSchema().find("ellipse.a").key
        keyB = refCat.getSchema().find("ellipse.b").key
        keyTheta = refCat.getSchema().find("ellipse.theta").key
        keyMag = refCat.getSchema().find("mag.%s" % exposure.getFilter().getName()).key
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
            outRecord.setPointD(self.keys["ref.center"], exposure.getWcs().skyToPixel(refRecord.getCoord()))

            # Next we determine the pixel region we want to fit.
            outRecord.setFootprint(setupFitRegion(self.config.fitRegion, exposure, srcRecord))

            # This is the WCS the parameters are defined in: it's a local tangent plane at
            # the position of the object, aligned with the celestial coordinate axes.
            # We have to use this coordinate system even when fitting to a single image
            # so the prior doesn't have to be transformed to the local coordinate system.
            # Note that the origin of this coordinate system is at the input position of
            # the object, so centroid parameters will be the *offset* from that position.
            fitWcs = self.config.makeFitWcs(outRecord.getCoord())

            # Now we'll transform the refCat ellipse to the fit coordinate system.  Because
            # the fit coordinate system is aligned with the coordinate axes, this should just
            # be a scaling and translation.
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

            # We now transform this ellipse and the refCat fluxes into the parameters defined by the model.
            ellipses = lsst.meas.multifit.Model.EllipseVector([ellipse2])
            self.model.readEllipses(
                ellipses,
                outRecord[self.keys["ref.nonlinear"]],
                outRecord[self.keys["ref.fixed"]]
                )
            # this flux->amplitudes conversion assumes the ref catalog is single-component, and that the
            # first component of the model is what that corresponds to; we may need to generalize this
            flux = self.fitCalib.getFlux(refRecord.get(keyMag))
            amplitudes = outRecord[self.keys["ref.amplitudes"]]
            amplitudes[:] = 0.0
            amplitudes[0] = flux

            # Finally, we tell the fitter to initialize the record, which includes setting "ref.parameters",
            # as only the fitter knows what parameters it will actually fit and how those map to
            # "ref.nonlinear" and "ref.amplitudes".  Depending on the fitter, it will probably attach
            # an initial PDF too.
            self.fitter.initialize(self.model, outRecord)

    def makeLikelihood(self, inputs, record):
        """Create a Likelihood object for a single object.

        The MeasureImage implementation creates a ProjectedLikelihood with data from a single
        exposure.
        """
        psf = TODO()
        return multifitLib.ProjectedLikelihood(
            self.model, record[self.keys["ref.fixed"]],
            self.config.makeFitWcs(record.getCoord()),
            self.fitCalib,
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
