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

import lsst.pex.config
import lsst.meas.base
from . import modelfitLib


class ShapeletPsfApproxConfig(lsst.pex.config.Config):
    models = lsst.pex.config.ConfigDictField(
        keytype=str,
        itemtype=modelfitLib.PsfFitterConfig,
        doc="a dictionary of models that can be used to fit the PSF",
        default={} # populated in setDefaults; can't do it on a single line here
        )
    sequence = lsst.pex.config.ListField(
        dtype=str,
        doc="a sequence of model names indicating which models should be fit, and their order",
        default=["DoubleShapelet"]
        )

    def setDefaults(self):
        super(ShapeletPsfApproxConfig, self).setDefaults()
        self.models["SingleGaussian"] = modelfitLib.PsfFitterConfig()
        self.models["SingleGaussian"].inner.order = -1
        self.models["SingleGaussian"].primary.order = 0
        self.models["SingleGaussian"].wings.order = -1
        self.models["SingleGaussian"].outer.order = -1
        self.models["DoubleGaussian"] = modelfitLib.PsfFitterConfig()
        self.models["DoubleGaussian"].inner.order = -1
        self.models["DoubleGaussian"].primary.order = 0
        self.models["DoubleGaussian"].wings.order = 0
        self.models["DoubleGaussian"].outer.order = -1
        self.models["DoubleShapelet"] = modelfitLib.PsfFitterConfig()
        self.models["DoubleShapelet"].inner.order = -1
        self.models["DoubleShapelet"].primary.order = 2
        self.models["DoubleShapelet"].wings.order = 1
        self.models["DoubleShapelet"].outer.order = -1
        self.models["Full"] = modelfitLib.PsfFitterConfig()
        self.models["Full"].inner.order = 0
        self.models["Full"].primary.order = 4
        self.models["Full"].wings.order = 4
        self.models["Full"].outer.order = 0

    def validate(self):
        super(ShapeletPsfApproxConfig, self).validate()
        if len(self.sequence) < 1:
            raise ValueError("sequence must have at least one element")
        for m in self.sequence:
            if m not in self.models:
                raise KeyError("All elements in sequence must be keys in models dict")

class ShapeletPsfApproxMixin(object):
    """Mixin base class for fitting shapelet approximations to the PSF model

    This class does almost all of the work for its two derived classes, ShapeletPsfApproxSingleFramePlugin
    and ShapeletPsfApproxForcedPlugin, which simply adapt it to the slightly different interfaces for
    single-frame and forced measurement.  It in turn delegates its work to the C++ PsfFitter class;
    it holds sequence of these corresponding to different models (generally with increasing complexity).
    Each PsfFitter starts with the result of the previous one as an input, using PsfFitter::adapt to
    hopefully allow these previous fits to reduce the time spent on the next one.
    """

    def __init__(self, config, name, schema):
        """Initialize the plugin, creating a sequence of PsfFitter instances to do the fitting and
        MultiShapeletFunctionKey instances to save the results to a record.
        """
        self.sequence = []
        for m in config.sequence:
            fitter = modelfitLib.PsfFitter(config.models[m].makeControl())
            modelKey = fitter.addModelFields(schema, schema.join(name, m))
            flagKey = schema.addField(schema.join(name, m, "flag"), type="Flag",
                                      doc="General failure flag set if anything goes wrong.")
            nIterKey = schema.addField(schema.join(name, m, "nIterations"), type=int,
                                       doc="Number of E-M steps taken by the fitter.")
            self.sequence.append((fitter, name, modelKey, flagKey, nIterKey))

    def measure(self, measRecord, exposure):
        """Fit the configured sequence of models the given Exposure's Psf, as evaluated at
        measRecord.getCentroid(), then save the results to measRecord.
        """
        if not exposure.hasPsf():
            raise lsst.meas.base.FatalAlgorithmError("ShapeletPsfApprox requires Exposure to have a Psf")
        psf = exposure.getPsf()
        psfImage = psf.computeKernelImage(measRecord.getCentroid())
        psfShape = psf.computeShape(measRecord.getCentroid())
        model = None
        lastFitter = None
        # Fit the first element in the sequence, using the PSFs moments to initialize the parameters
        # For every other element in the fitting sequence, use the previous fit to initialize the parameters
        for fitter, name, modelKey, flagKey, nIterKey in self.sequence:
            try:
                if model is None:
                    model = fitter.makeInitial(psfShape)
                else:
                    model = fitter.adapt(model, lastFitter)
                lastFitter = fitter
                nIterations = fitter.apply(model, psfImage)
                measRecord.set(modelKey, model)
                measRecord.set(nIterKey, nIterations)
                if nIterations == fitter.getMaxIterations():
                    measRecord.set(flagKey, True)
            except Exception:
                # Since this is an unexpected failure mode, we throw now and do not proceed to the next fit.
                measRecord.set(flagKey, True)
                raise

    # This plugin doesn't need to set a flag in fail(), because it should have been
    # done already by the try/except in measure().
    def fail(self, measRecord, error=None):
        pass

class ShapeletPsfApproxSingleFrameConfig(lsst.meas.base.SingleFramePluginConfig, ShapeletPsfApproxConfig):

    def setDefaults(self):
        lsst.meas.base.SingleFramePluginConfig.setDefaults(self)
        ShapeletPsfApproxConfig.setDefaults(self)

@lsst.meas.base.register("modelfit_ShapeletPsfApprox")
class ShapeletPsfApproxSingleFramePlugin(lsst.meas.base.SingleFramePlugin, ShapeletPsfApproxMixin):
    """Minimal subclass of ShapeletPsfApproxMixin to conform to the single-frame measurement API.

    This class simply provides __init__ and measure methods that matched the SingleFramePlugin signatures
    and delegate to the ShapeletPsfApproxMixin's implmeentations.
    """
    ConfigClass = ShapeletPsfApproxSingleFrameConfig

    @staticmethod
    def getExecutionOrder():
        return 1.0

    def __init__(self, config, name, schema, metadata):
        ShapeletPsfApproxMixin.__init__(self, config, name, schema)
        lsst.meas.base.SingleFramePlugin.__init__(self, config, name, schema, metadata)

    def measure(self, measRecord, exposure):
        ShapeletPsfApproxMixin.measure(self, measRecord, exposure)

    def fail(self, measRecord, error=None):
        ShapeletPsfApproxMixin.fail(self, measRecord, error)

class ShapeletPsfApproxForcedConfig(lsst.meas.base.ForcedPluginConfig, ShapeletPsfApproxConfig):

    def setDefaults(self):
        lsst.meas.base.ForcedPluginConfig.setDefaults(self)
        ShapeletPsfApproxConfig.setDefaults(self)

@lsst.meas.base.register("modelfit_ShapeletPsfApprox")
class ShapeletPsfApproxForcedPlugin(lsst.meas.base.ForcedPlugin, ShapeletPsfApproxMixin):
    """Minimal subclass of ShapeletPsfApproxMixin to conform to the forced measurement API.

    This class simply provides __init__ and measure methods that matched the ForcedPlugin signatures
    and delegate to the ShapeletPsfApproxMixin's implmeentations.
    """
    ConfigClass = ShapeletPsfApproxForcedConfig

    @staticmethod
    def getExecutionOrder():
        return 1.0

    def __init__(self, config, name, schemaMapper, metadata):
        ShapeletPsfApproxMixin.__init__(self, config, name, schemaMapper.editOutputSchema())
        lsst.meas.base.ForcedPlugin.__init__(self, config, name, schemaMapper, metadata)

    def measure(self, measRecord, exposure, refRecord, refWcs):
        ShapeletPsfApproxMixin.measure(self, measRecord, exposure)

    def fail(self, measRecord, error=None):
        ShapeletPsfApproxMixin.fail(self, measRecord, error)
