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

__all__ = ("BaseMeasureConfig", "BaseMeasureTask")

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

    def setDefaults(self):
        self.psf.innerOrder = 4
        self.psf.outerOrder = 0

class BaseMeasureTask(lsst.pipe.base.CmdLineTask):
    """Base class for MeasureImageTask and MeasureMultiTask to aggregate shared code"""

    def addFields(self, prefix, doc):
        """Add <prefix>.ellipse and <prefix>.center keys to self.schema, saving keys in self.keys
        @param[in] prefix       key name prefix
        @param[in] doc          documentation prefix
        """
        self.keys["%s.ellipse" % prefix] = self.schema.addField("%s.ellipse" % prefix, type="MomentsD",
                                                                doc=("%s ellipse" % doc))
        self.keys["%s.center" % prefix] = self.schema.addField("%s.center" % prefix, type="PointD",
                                                               doc=("%s center position" % doc))
    def addDerivedFields(self):
        """Add fields to self.schema for quantities derived from the SampleSet
        """
        self.addFields("mean", "Posterior mean")
        self.addFields("median", "Posterior median")

    def fillDerivedFields(self, record):
        samples = record.getSamples()
        parameterDef = samples.getParameterDefinition()
        mean = parameterDef.makeEllipse(samples.computeMean(), record.getPointD(self.keys["source.center"]))
        record.set(self.keys["mean.ellipse"], lsst.afw.geom.ellipses.Quadrupole(mean.getCore()))
        record.set(self.keys["mean.center"], mean.getCenter())
        median = parameterDef.makeEllipse(samples.computeQuantiles(numpy.array([0.5])),
                                          record.getPointD(self.keys["source.center"]))
        record.set(self.keys["median.ellipse"], lsst.afw.geom.ellipses.Quadrupole(median.getCore()))
        record.set(self.keys["median.center"], median.getCenter())
