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
        dtype=multifitLib.UnitTransformedLikelihood.ConfigClass,
        doc="Config for likelihood object that computes model probability at given parameters"
    )
    minInitialRadius = lsst.pex.config.Field(
        dtype=float,
        default=0.1,
        doc="Minimum deconvolved initial radius in pixels"
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
        srcCat = inputs.srcCat

        exposureUnits = multifitLib.UnitSystem(inputs.exposure)
        exposurePsf = inputs.exposure.getPsf()
        exposureCalib = inputs.exposure.getCalib()

        # SchemaMapper will transfer ID, Coord, Footprint (we'll overwrite the Coord with that from RefCat)
        mapper = lsst.afw.table.SchemaMapper(lsst.afw.table.SourceTable.makeMinimalSchema())
        mapper.addMinimalSchema(lsst.meas.multifit.ModelFitTable.makeMinimalSchema())

        for srcRecord in srcCat:
            if (srcRecord.getShapeFlag() or srcRecord.getCentroidFlag()
                or srcRecord.getApFluxFlag() or srcRecord.getPsfFluxFlag()):
                continue

            outRecord = outCat.addNew()

            # Start by setting some miscellaneous calculated fields
            outRecord.assign(srcRecord, mapper)
            outRecord.setD(self.keys["snr"], srcRecord.getApFlux() / srcRecord.getApFluxErr())
            outRecord.setCoord(srcRecord.getCoord())
            outRecord.setPointD(self.keys["center"], srcRecord.getCentroid())

            # Next we determine the pixel region we want to fit.
            outRecord.setFootprint(setupFitRegion(self.config.fitRegion, inputs.exposure, srcRecord))

            # Setup the coordinate and photometric systems to use for the parameters, and the transform
            # from the exposure to the parameter unit system
            nominalMag = exposureCalib.getMagnitude(srcRecord.getPsfFlux())
            units = self.makeUnitSystem(outRecord, outRecord.getCoord(), nominalMag)
            transform = multifitLib.LocalUnitTransform(outRecord.getCoord(), exposureUnits, units)

            # Start with the ellipse from the Shape and Centroid src slots (should refine this),
            # subtract the PSF moments (truncate radius as specified by config), and transform
            # to the parameter unit system
            fullShape = srcRecord.getShape()
            psfShape = exposurePsf.computeShape(srcRecord.getCentroid())
            deconvolvedShape = lsst.afw.geom.ellipses.Quadrupole(
                max(fullShape.getIxx() - psfShape.getIxx(), self.config.minInitialRadius**2),
                max(fullShape.getIyy() - psfShape.getIyy(), self.config.minInitialRadius**2),
                fullShape.getIxy() - psfShape.getIxy()
                )
            ellipse1 = lsst.afw.geom.ellipses.Ellipse(deconvolvedShape, srcRecord.getCentroid())
            ellipse2 = ellipse1.transform(transform.geometric)

            # Fill the initial.nonlinear and fixed parameter arrays from the initial ellipse
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

        The MeasureImage implementation creates a UnitTransformedLikelihood with data from a single
        exposure.
        """
        psfModel = FitPsfAlgorithm.apply(self.config.psf.makeControl(), inputs.exposure.getPsf(),
                                         record.get(self.keys["center"]))
        psf = psfModel.asMultiShapelet()
        return multifitLib.UnitTransformedLikelihood(
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
