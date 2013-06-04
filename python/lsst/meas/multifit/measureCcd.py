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

from .measureImage import MeasureImageTask

__all__ = ("MeasureCcdConfig", "MeasureCcdTask")

MeasureCcdConfig = MeasureImageTask.ConfigClass

class MeasureCcdTask(MeasureImageTask):
    """Specialization of MeasureImageTask for running on calexps, after processCcd.py or processEimage.py
    """

    _DefaultName = "measureCcd"

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser
        """
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "calexp", help="data ID, e.g. --id visit=1 raft=2,2 sensor=1,1")
        return parser
