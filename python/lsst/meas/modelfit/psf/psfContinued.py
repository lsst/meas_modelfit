#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2014 LSST Corporation.
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

from __future__ import absolute_import, division, print_function

# all new classes here are accessed via registries, not direct imports.
__all__ = (
    "GeneralPsfFitterComponentConfig",
    "GeneralPsfFitterConfig"
)

import lsst.pex.config
import lsst.meas.base
from .psf import (
    GeneralPsfFitterControl, GeneralPsfFitterComponentControl,
    GeneralPsfFitter, GeneralPsfFitterAlgorithm,
    DoubleShapeletPsfApproxAlgorithm, DoubleShapeletPsfApproxControl
)


lsst.meas.base.wrapSimpleAlgorithm(
    DoubleShapeletPsfApproxAlgorithm,
    Control=DoubleShapeletPsfApproxControl,
    module='lsst.meas.modelfit',
    name='modelfit_DoubleShapeletPsfApprox',
    executionOrder=lsst.meas.base.BasePlugin.SHAPE_ORDER
)


GeneralPsfFitterComponentConfig = lsst.pex.config.makeConfigClass(
    GeneralPsfFitterComponentControl,
    module='lsst.meas.modelfit'
)
GeneralPsfFitterConfig = lsst.pex.config.makeConfigClass(
    GeneralPsfFitterControl,
    module='lsst.meas.modelfit'
)
GeneralPsfFitter.ConfigClass = GeneralPsfFitterConfig


class GeneralShapeletPsfApproxConfig(lsst.pex.config.Config):
    models = lsst.pex.config.ConfigDictField(
        keytype=str,
        itemtype=GeneralPsfFitterConfig,
        doc="a dictionary of models that can be used to fit the PSF",
        default={} # populated in setDefaults; can't do it on a single line
    )
    sequence = lsst.pex.config.ListField(
        dtype=str,
        doc=("a sequence of model names indicating which models should be fit,"
             " and their order"),
        default=["DoubleShapelet"]
    )

    def setDefaults(self):
        super(GeneralShapeletPsfApproxConfig, self).setDefaults()
        self.models["SingleGaussian"] = GeneralPsfFitterConfig()
        self.models["SingleGaussian"].inner.order = -1
        self.models["SingleGaussian"].primary.order = 0
        self.models["SingleGaussian"].wings.order = -1
        self.models["SingleGaussian"].outer.order = -1
        self.models["DoubleGaussian"] = GeneralPsfFitterConfig()
        self.models["DoubleGaussian"].inner.order = -1
        self.models["DoubleGaussian"].primary.order = 0
        self.models["DoubleGaussian"].wings.order = 0
        self.models["DoubleGaussian"].outer.order = -1
        self.models["DoubleShapelet"] = GeneralPsfFitterConfig()
        self.models["DoubleShapelet"].inner.order = -1
        self.models["DoubleShapelet"].primary.order = 2
        self.models["DoubleShapelet"].wings.order = 1
        self.models["DoubleShapelet"].outer.order = -1
        self.models["Full"] = GeneralPsfFitterConfig()
        self.models["Full"].inner.order = 0
        self.models["Full"].primary.order = 4
        self.models["Full"].wings.order = 4
        self.models["Full"].outer.order = 0

    def validate(self):
        super(GeneralShapeletPsfApproxConfig, self).validate()
        if len(self.sequence) < 1:
            raise ValueError("sequence must have at least one element")
        for m in self.sequence:
            if m not in self.models:
                raise KeyError(
                    "All elements in sequence must be keys in models dict"
                )


class GeneralShapeletPsfApproxMixin(object):
    """Mixin base class for fitting shapelet approximations to the PSF model

    This class does almost all of the work for its two derived classes,
    GeneralShapeletPsfApproxSingleFramePlugin and
    GeneralShapeletPsfApproxForcedPlugin, which simply adapt it to the
    slightly different interfaces for single-frame and forced measurement.  It
    in turn delegates its work to the C++ GeneralPsfFitter class; it holds
    sequence of these corresponding to different models (generally with
    increasing complexity). Each GeneralPsfFitter starts with the result of
    the previous one as an input, using GeneralPsfFitter::adapt to hopefully
    allow these previous fits to reduce the time spent on the next one.

    At present, this plugin does not define any failure flags, which will
    almost certainly have to be changed in the future.  So far, however, I
    haven't actually seen it fail on any PSFs I've given it, so I'll wait
    until we can run on large enough data volumes to see what the actual
    failure modes are, instead of trying to guess them in advance.
    """

    def __init__(self, config, name, schema):
        """Initialize the plugin, creating a sequence of GeneralPsfFitter
        instances to do the fitting and MultiShapeletFunctionKey instances to
        save the results to a record.
        """
        self.sequence = []
        for m in config.sequence:
            fitter = GeneralPsfFitterAlgorithm(
                config.models[m].makeControl(),
                schema,
                schema[name][m].getPrefix()
            )
            self.sequence.append((fitter, schema[name][m].getPrefix()))

    def measure(self, measRecord, exposure):
        """Fit the configured sequence of models the given Exposure's Psf, as
        evaluated at measRecord.getCentroid(), then save the results to
        measRecord.
        """
        if not exposure.hasPsf():
            raise lsst.meas.base.FatalAlgorithmError(
                "GeneralShapeletPsfApprox requires Exposure to have a Psf")
        psf = exposure.getPsf()
        psfImage = psf.computeKernelImage(measRecord.getCentroid())
        psfShape = psf.computeShape(measRecord.getCentroid())
        lastError = None
        lastModel = None
        # Fit the first element in the sequence, using the PSFs moments to
        # initialize the parameters For every other element in the fitting
        # sequence, use the previous fit to initialize the parameters
        for fitter, name in self.sequence:
            try:
                if lastModel is None:
                    fitter.measure(measRecord, psfImage, psfShape)
                else:
                    fitter.measure(measRecord, psfImage,
                                   fitter.adapt(lastResult, lastModel))
                lastResult = measRecord.get(fitter.getKey())
                lastModel = fitter.getModel()
            except lsst.meas.base.baseMeasurement.FATAL_EXCEPTIONS:
                raise
            except lsst.meas.base.MeasurementError as error:
                fitter.fail(measRecord, error.cpp)
                lastError = error
            except Exception as error:
                fitter.fail(measRecord)
                lastError = error
        # When we are done with all the fitters, raise the last error if there
        # was one. This gives the calling task a chance to do whatever it
        # wants
        if not lastError is None:
            raise lastError

    # This plugin doesn't need to set a flag on fail, because it should have
    # been done already by the individual fitters in the sequence
    def fail(self, measRecord, error=None):
        pass


class GeneralShapeletPsfApproxSingleFrameConfig(
    lsst.meas.base.SingleFramePluginConfig,
    GeneralShapeletPsfApproxConfig
):

    def setDefaults(self):
        lsst.meas.base.SingleFramePluginConfig.setDefaults(self)
        GeneralShapeletPsfApproxConfig.setDefaults(self)


@lsst.meas.base.register("modelfit_GeneralShapeletPsfApprox")
class GeneralShapeletPsfApproxSingleFramePlugin(
    lsst.meas.base.SingleFramePlugin,
    GeneralShapeletPsfApproxMixin
):
    """Minimal subclass of GeneralShapeletPsfApproxMixin to conform to the
    single-frame measurement API.

    This class simply provides __init__ and measure methods that matched the
    SingleFramePlugin signatures and delegate to the
    GeneralShapeletPsfApproxMixin's implementations.
    """
    ConfigClass = GeneralShapeletPsfApproxSingleFrameConfig

    @staticmethod
    def getExecutionOrder():
        return 1.0

    def __init__(self, config, name, schema, metadata):
        GeneralShapeletPsfApproxMixin.__init__(self, config, name, schema)
        lsst.meas.base.SingleFramePlugin.__init__(self, config, name, schema,
                                                  metadata)

    def measure(self, measRecord, exposure):
        GeneralShapeletPsfApproxMixin.measure(self, measRecord, exposure)

    def fail(self, measRecord, error=None):
        GeneralShapeletPsfApproxMixin.fail(self, measRecord, error)


class GeneralShapeletPsfApproxForcedConfig(
    lsst.meas.base.ForcedPluginConfig,
    GeneralShapeletPsfApproxConfig
):

    def setDefaults(self):
        lsst.meas.base.ForcedPluginConfig.setDefaults(self)
        GeneralShapeletPsfApproxConfig.setDefaults(self)


@lsst.meas.base.register("modelfit_GeneralShapeletPsfApprox")
class GeneralShapeletPsfApproxForcedPlugin(
    lsst.meas.base.ForcedPlugin,
    GeneralShapeletPsfApproxMixin
):
    """Minimal subclass of GeneralShapeletPsfApproxMixin to conform to the
    forced measurement API.

    This class simply provides __init__ and measure methods that matched the
    ForcedPlugin signatures and delegate to the
    GeneralShapeletPsfApproxMixin's implementations.
    """
    ConfigClass = GeneralShapeletPsfApproxForcedConfig

    @staticmethod
    def getExecutionOrder():
        return 1.0

    def __init__(self, config, name, schemaMapper, metadata):
        GeneralShapeletPsfApproxMixin.__init__(self, config, name,
                                               schemaMapper.editOutputSchema())
        lsst.meas.base.ForcedPlugin.__init__(self, config, name, schemaMapper,
                                             metadata)

    def measure(self, measRecord, exposure, refRecord, refWcs):
        GeneralShapeletPsfApproxMixin.measure(self, measRecord, exposure)

    def fail(self, measRecord, error=None):
        GeneralShapeletPsfApproxMixin.fail(self, measRecord, error)
