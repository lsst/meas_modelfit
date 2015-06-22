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
from . import modelfitLib

class OptFitStageConfig(lsst.pex.config.Config):
    model = modelRegistry.makeField(doc="Model to use for this stage.", default="bulge+disk")
    prior = priorRegistry.makeField(doc="Prior to use for this stage.", default="linear")
    optimizer = lsst.pex.config.ConfigField(dtype=modelfitLib.OptimizerConfig,
                                            doc="Optimizer configuration for this stage.")
    psf = lsst.pex.config.Field(dtype=str, default="full",
                                doc="Name of the ShapeletPsfApprox model to use in this stage.")

class OptFitStage(object):

    def __init__(self, config, algName, stageName, schema):
        self.config = config
        self.model = self.config.model.apply()
        self.prior = self.config.prior.apply()
        self.keys = {}
        self.keys["psf"] = lsst.shapelet.MultiShapeletFunctionKey(
            schema["modelfit"]["ShapeletPsfApprox"][config.psf]
        )
        def addField(name, type, doc):
            return schema.addField(
                "{algName}_{stageName}_{name}".format(algName=algName, stageName=stageName, name=name),
                type=type,
                doc=doc
            )
        self.keys["nonlinear"] = lsst.afw.table.ArrayDKey([
            addField(name, numpy.float64, "best-fit nonlinear parameter {}".format(name))
            for name in self.model.getNonlinearNames()
        ])
        self.keys["amplitudes"] = lsst.afw.table.ArrayDKey([
            addField(name, numpy.float64, "best-fit amplitude parameter {}".format(name))
            for name in self.model.getAmplitudeNames()
        ])
        self.keys["fixed"] = lsst.afw.table.ArrayDKey([
            addField(name, numpy.float64, "fixed parameter {} (not fit)".format(name))
            for name in self.model.getFixedNames()
        ])
        self.keys["flag"] = {}
        self.keys["flag"]["failure"] = addField(
            "flag", "Flag", "general failure flag set if model fit failed"
        )
        self.keys["flag"]["trSmall"] = addField(
            "flag_trSmall", "Flag", "trust region grew too small before gradient threshold was met"
        )
        self.keys["flag"]["maxIter"] = addField(
            "flag_maxIter", "Flag", "fit exceeded the maximum number of allowed iterations"
        )
        self.keys["flag"]["numericError"] = addField(
            "flag_numericError", "Flag", "numeric error (usually overflow, underflow, or NaNs)"
        )


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
        default=["initial", "final"]
    )

    def setDefaults(self):
        lsst.meas.base.SingleFramePluginConfig.setDefaults(self)
        self.stages["initial"] = OptFitStageConfig()
        self.stages["initial"].model = "gaussian"
        self.stages["initial"].psf = "DoubleGaussian"
        self.stages["initial"].optimizer.minTrustRadiusThreshold = 1E-2
        self.stages["initial"].optimizer.gradientThreshold = 1E-2
        self.stages["final"] = OptFitStageConfig()
        self.stages["final"].model = "bulge+disk"
        self.stages["final"].psf = "Full"


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
