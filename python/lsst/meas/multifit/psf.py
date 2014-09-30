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

class PsfFitterBasePluginConfig(lsst.meas.base.BasePluginConfig):
    models = lsst.pex.config.ConfigDictDict(
        dtype=multifitLib.PsfFitterConfig,
        doc="a dictionary of models that can be used to fit the PSF",
        default={
            "SingleGaussian": multifitLib.PsfFitterConfig()
            "DoubleGaussian": multifitLib.PsfFitterConfig(),
            "Full": multifitLib.PsfFitterConfig()
            }
        )
    sequence = lsst.pex.config.ListField(
        dtype=str,
        doc="a sequence of model names indicating which models should be fit, and their order",
        default=["DoubleGaussian"]
        )

    def setDefaults(self):
        super(PsfFitterPluginConfig, self).setDefaults()
        self.models["SingleGaussian"].inner.order = -1
        self.models["SingleGaussian"].primary.order = 0
        self.models["SingleGaussian"].wings.order = -1
        self.models["SingleGaussian"].outer.order = -1
        self.models["DoubleGaussian"].inner.order = -1
        self.models["DoubleGaussian"].primary.order = 0
        self.models["DoubleGaussian"].wings.order = 0
        self.models["DoubleGaussian"].outer.order = -1
        self.models["Full"].inner.order = 0
        self.models["Full"].primary.order = 4
        self.models["Full"].wings.order = 4
        self.models["Full"].outer.order = 0

class PsfFitterSingleFramePluginConfig(lsst.meas.base.SingleFramePluginConfig, PsfFitterBasePluginConfig):
    pass
