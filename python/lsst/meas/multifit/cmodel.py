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
from . import multifitLib

class CModelSingleFrameConfig(lsst.meas.base.SingleFramePluginConfig, multifitLib.CModelConfig):

    def setDefaults(self):
        lsst.meas.base.SingleFramePluginConfig.setDefaults(self)
        multifitLib.CModelConfig.setDefaults(self)
        self.executionOrder = 3.0


@lsst.meas.base.register("multifit_CModel")
class CModelSingleFramePlugin(lsst.meas.base.SingleFramePlugin):
    """Single-frame measurement interface for CModelAlgorithm.

    This class simply provides __init__ and measure methods that matched the SingleFramePlugin signatures
    and delegate to the CModelAlgorithm's methods.
    """
    ConfigClass = CModelSingleFrameConfig

    def __init__(self, config, name, schema, flags, others, metadata):
        lsst.meas.base.SingleFramePlugin.__init__(self, config, name, schema, flags, others, metadata)
        self.algorithm = multifitLib.CModelAlgorithm(name, config.makeControl(), schema)

    def measure(self, measRecord, exposure):
        self.algorithm.measure(measRecord, exposure)

    def fail(self, measRecord, error=None):
        self.algorithm.fail(measRecord, error.cpp if error is not None else None)

class CModelForcedConfig(lsst.meas.base.ForcedPluginConfig, multifitLib.CModelConfig):

    def setDefaults(self):
        lsst.meas.base.ForcedPluginConfig.setDefaults(self)
        multifitLib.CModelConfig.setDefaults(self)
        self.executionOrder = 3.0

@lsst.meas.base.register("multifit_CModel")
class CModelForcedPlugin(lsst.meas.base.ForcedPlugin):
    """Minimal subclass of CModelMixin to conform to the forced measurement API.

    This class simply provides __init__ and measure methods that matched the ForcedPlugin signatures
    and delegate to the CModelMixin's implmentations.
    """
    ConfigClass = CModelForcedConfig

    def __init__(self, config, name, schemaMapper, flags, others, metadata):
        lsst.meas.base.ForcedPlugin.__init__(self, config, name, schemaMapper, flags, others, metadata)
        self.algorithm = multifitLib.CModelAlgorithm(name, config.makeControl(), schemaMapper)

    def measure(self, measRecord, exposure, refRecord, refWcs):
        self.algorithm.measure(measRecord, exposure)

    def fail(self, measRecord, error=None):
        self.algorithm.fail(measRecord, error.cpp if error is not None else None)
