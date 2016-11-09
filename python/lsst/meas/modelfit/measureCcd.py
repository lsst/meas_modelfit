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

import lsst.pipe.base
import lsst.pex.config

from .measureImage import MeasureImageTask

try:
    from lsst.meas.mosaic import applyMosaicResults
except ImportError:
    applyMosaicResults = None

__all__ = ("MeasureCcdConfig", "MeasureCcdTask")


class MeasureCcdConfig(MeasureImageTask.ConfigClass):
    doApplyUberCal = lsst.pex.config.Field(
        dtype=bool,
        doc="Apply meas_mosaic ubercal results to input calexps?",
        default=False
    )


class MeasureCcdTask(MeasureImageTask):
    """Specialization of MeasureImageTask for running on calexps, after processCcd.py or processEimage.py
    """

    ConfigClass = MeasureCcdConfig

    _DefaultName = "measureCcd"

    def readInputs(self, dataRef):
        inputs = MeasureImageTask.readInputs(self, dataRef)
        if self.config.doApplyUberCal:
            if not applyMosaicResults:
                raise RuntimeError(
                    "Cannot use improved calibrations for %s because meas_mosaic could not be imported."
                    % dataRef.dataId
                )
            applyMosaicResults(dataRef, calexp=inputs.exposure)
        return inputs

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser
        """
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "calexp", help="data ID, e.g. --id visit=1 raft=2,2 sensor=1,1")
        return parser

    def getPreviousTaskClass(self):
        return MeasureCcdTask
