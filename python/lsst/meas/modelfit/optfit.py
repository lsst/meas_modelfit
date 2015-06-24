#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2015 LSST/AURA.
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
import lsst.meas.base
import lsst.afw.table
import lsst.shapelet
import lsst.afw.geom.ellipses

from .models import modelRegistry
from .priors import priorRegistry
from . import fitRegion
from . import modelfitLib

class OptFitStageConfig(lsst.pex.config.Config):
    model = modelRegistry.makeField(doc="Model to use for this stage.", default="bulge+disk")
    prior = priorRegistry.makeField(doc="Prior to use for this stage.", default="linear")
    optimizer = lsst.pex.config.ConfigField(dtype=modelfitLib.OptimizerConfig,
                                            doc="Optimizer configuration for this stage.")
    region = lsst.pex.config.ConfigField(dtype=fitRegion.FitRegionConfig,
                                            doc="Configuration for which pixels to include in the fit.")
    psf = lsst.pex.config.Field(dtype=str, default="Full",
                                doc="Name of the ShapeletPsfApprox model to use in this stage.")

class OptFitStage(object):

    def __init__(self, config, algName, stageName, schema):
        self.config = config
        self.model = self.config.model.apply()
        self.prior = self.config.prior.apply()
        self.keys = {}
        # Add keys for the shapelet PSF approximation that's computed by another plugin.
        self.keys["psf"] = lsst.shapelet.MultiShapeletFunctionKey(
            schema["modelfit"]["ShapeletPsfApprox"][config.psf]
        )
        # Some local convenience functions to add fields for outputs and save the keys
        def addArrayFields(aggregate, names, doc):
            self.keys[aggregate] = lsst.afw.table.ArrayDKey([
                schema.addField("_".join([algName, stageName, name]), type=numpy.float64,
                                doc=doc.format(name))
                for name in names
            ])
        def addFlagField(name, doc):
            f = [algName, stageName, "flag"]
            if name is not None: f.append(name)
            self.keys["flag"][name] = schema.addField("_".join(f), type="Flag", doc=doc)
        # Add fields for outputs, and save their keys
        addArrayFields("nonlinear", self.model.getNonlinearNames(), "best-fit nonlinear parameter {}")
        addArrayFields("amplitudes", self.model.getAmplitudeNames(), "best-fit amplitude parameter {}")
        addArrayFields("fixed", self.model.getFixedNames(), "fixed parameter {} (not fit)")
        self.keys["flag"] = {}
        addFlagField(None, "general failure flag set if model fit failed")
        addFlagField("trSmall", "trust region grew too small before gradient threshold was met")
        addFlagField("maxIter", "fit exceeded the maximum number of allowed iterations")
        addFlagField("numericError", "numeric error (usually overflow, underflow, or NaNs)")


class OptFitConfig(lsst.meas.base.SingleFramePluginConfig):
    stages = lsst.pex.config.ConfigDictField(
        keytype=str,
        itemtype=OptFitStageConfig,
        doc="a dictionary of stages of model fitting that can be executed in sequence",
        default={} # populated in setDefaults; can't do it on a single line here
    )
    sequence = lsst.pex.config.ListField(
        dtype=str,
        doc="a sequence of stage names indicating which models should be fit, and their order",
        default=["initial", "tiedBD"]
    )


    def setDefaults(self):
        lsst.meas.base.SingleFramePluginConfig.setDefaults(self)
        # A fast initial stage that just fits a Gaussian model with a DoubleGaussian PSF approximation;
        # should be used to warm-start other slower models.
        self.stages["initial"] = OptFitStageConfig()
        self.stages["initial"].model.name = "gaussian"
        self.stages["initial"].psf = "DoubleGaussian"
        self.stages["initial"].optimizer.minTrustRadiusThreshold = 1E-2
        self.stages["initial"].optimizer.gradientThreshold = 1E-2
        # A bulge+disk model where the two components have the same ellipticity and a fixed radius ratio
        # (as used in e.g. lensfit).
        self.stages["tiedBD"] = OptFitStageConfig()
        self.stages["tiedBD"].model.name = "bulge+disk"
        # --------------------------------------------------------------------------------------------
        # The stages below aren't enabled by default (they're not in the "sequence" field's defaults),
        # but they are sufficiently useful we want to make them easy for users to enable.
        # --------------------------------------------------------------------------------------------
        # An exponential model, using the SDSS approximation.
        self.stages["exp"] = OptFitStageConfig()
        self.stages["exp"].model.name = "fixed-sersic"
        self.stages["exp"].model["fixed-sersic"].profile = "lux"
        # A de Vaucouleur model, using the SDSS approximation.
        self.stages["dev"] = OptFitStageConfig()
        self.stages["dev"].model.name = "fixed-sersic"
        self.stages["dev"].model["fixed-sersic"].profile = "luv"
        # A full bulge+disk model with independent ellipses (but still shared center)
        self.stages["fullBD"] = OptFitStageConfig()
        self.stages["fullBD"].model.name = "bulge+disk"
        self.stages["fullBD"].model["bulge+disk"].bulgeRadius = None


@lsst.meas.base.register("modelfit_OptFit")
class OptFitPlugin(lsst.meas.base.SingleFramePlugin):
    """Single-frame measurement plugin for greedy fitting of models to sources.
    """
    ConfigClass = OptFitConfig

    @staticmethod
    def getExecutionOrder():
        return 3.0

    def __init__(self, config, name, schema, metadata):
        lsst.meas.base.SingleFramePlugin.__init__(self, config, name, schema, metadata)
        self.stages = [OptFitStage(self.config.stages[stageName], name, stageName, schema)
                       for stageName in self.config.sequence]

    def measure(self, measRecord, exposure):
        pass

    def fail(self, measRecord, error=None):
        self.algorithm.fail(measRecord, error.cpp if error is not None else None)
