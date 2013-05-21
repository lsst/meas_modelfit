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

import lsst.pex.config
import lsst.pipe.base
import lsst.afw.table
import lsst.meas.extensions.multiShapelet

from .samplers import BaseSamplerTask, NaiveGridSamplerTask
from .multifitLib import SingleEpochObjective, ModelFitCatalog, ModelFitTable
from .models import BulgeDiskModelConfig
from .fitRegion import setupFitRegion

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
    model = lsst.pex.config.ConfigurableField(
        target=BulgeDiskModelConfig.makeBasis,
        doc="Definition of the galaxy model to fit"
    )
    psf = lsst.pex.config.ConfigField(
        dtype=lsst.meas.extensions.multiShapelet.FitPsfConfig,
        doc="Config options for approximating the PSF using shapelets"
    )
    fitRegion = lsst.pex.config.ConfigField(
        dtype=setupFitRegion.ConfigClass,
        doc="Parameters that control which pixels to include in the model fit"
    )

class MeasureImageTask(lsst.pipe.base.CmdLineTask):
    ConfigClass = MeasureImageConfig
    dataPrefix = ""

    def __init__(self, **kwds):
        lsst.pipe.base.CmdLineTask.__init__(self, **kwds)
        self.schemaMapper = lsst.afw.table.SchemaMapper(lsst.afw.table.SourceTable.makeMinimalSchema())
        self.schemaMapper.addMinimalSchema(lsst.meas.multifit.ModelFitTable.makeMinimalSchema())
        self.schema = self.schemaMapper.getOutputSchema()
        self.fitPsf = self.config.psf.makeControl().makeAlgorithm(self.schema)
        self.makeSubtask("sampler", schema=self.schema, model=self.config.model)
        self.basis = self.model.apply()
        self.ObjectiveClass = SingleEpochObjective[self.basis.getSize()]

    def readInputs(self, dataRef):
        """Return a lsst.pipe.base.Struct containing:
          - exposure ----- lsst.afw.image.ExposureF to fit
          - catalog ------ lsst.afw.table.SourceCatalog with initial measurements
        """
        return lsst.pipe.base.Struct(
            exposure = dataRef.get(self.dataPrefix + "calexp", immediate=True),
            catalog = dataRef.get(self.dataPrefix + "src", immediate=True),
            )

    def writeOutputs(self, dataRef, outputs):
        """Write task outputs using the butler.
        """
        dataRef.put(outputs.catalog, self.dataPrefix + "modelfits")

    def processObject(self, exposure, source, record):
        """Process a single object.

        @param[in] exposure    lsst.afw.image.ExposureF to fit
        @param[in] source      lsst.afw.table.SourceRecord that defines the object to be measured.
                               Must have valid "shape" and "centroid" slots.
        @param[in] record      ModelFitRecord with fields for output.
        """
        if not exposure.hasPsf():
            raise RuntimeError("Exposure has no PSF")
        psfModel = self.fitPsf.apply(record, exposure.getPsf(), source.getCentroid())
        if psfModel.hasFailed():
            raise RuntimeError("Shapelet approximation to PSF failed")
        psf = psfModel.asMultiShapelet()
        footprint = setupFitRegion(self.config.fitRegion, exposure, source)
        sampler = self.sampler.setup(exposure=inputs.exposure, source=source)
        objective = self.ObjectiveClass(
            self.config.objective.makeControl(), basis, psf,
            exposure.getMaskedImage(), footprint
        )
        samples = sampler.run(objective)
        # TODO: apply prior here!
        record.setSamples(samples)
        ellipse = sampler.interpret(samples.computeMean())
        record.setCentroid(ellipse.getCenter())
        record.setShape(ellipse.getCore())

    def processImage(self, exposure, catalog):
        """Process all sources in an exposure, drawing samples from each sources, returning
        the outputs as a Struct with a single 'catalog' attribute.
        """
        outputs = lsst.pipe.base.Struct(
            catalog = ModelFitCatalog(self.schema)
        )
        for source in inputs.catalog:
            record = outputs.catalog.addNew()
            record.assign(source, self.schemaMapper)
            self.processObject(exposure=inputs.exposure, source=source, record=record)
            samples.write(results)
        return outputs

    def run(self, dataRef):
        """Process the exposure/catalog associated with the given dataRef, and write the
        output catalog using the butler.
        """
        inputs = self.readInputs(dataRef)
        outputs = self.processImage(exposure=inputs.exposure, catalog=inputs.catalog)
        self.writeOutputs(dataRef, outputs)
        return outputs

    def getSchemaCatalogs(self):
        """Return a dict of empty catalogs for each catalog dataset produced by this task."""
        return {self.dataPrefix + "modelfits": ModelFitCatalog(self.schema)}
