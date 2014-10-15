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
from lsst.pipe.tasks.coaddBase import CoaddDataIdContainer

from .measureImage import MeasureImageTask

__all__ = ("MeasureCoaddConfig", "MeasureCoaddTask")

class MeasureCoaddConfig(MeasureImageTask.ConfigClass):
    """Config for ProcessCoadd"""
    coaddName = lsst.pex.config.Field(
        doc = "coadd name: typically one of deep or goodSeeing",
        dtype = str,
        default = "deep",
    )

class MeasureCoaddTask(MeasureImageTask):
    """Specialization of MeasureImageTask for running on coadds, after processCoadd.py
    """

    _DefaultName = "measureCoadd"
    ConfigClass = MeasureCoaddConfig

    def __init__(self, **kwargs):
        MeasureImageTask.__init__(self, **kwargs)
        self.dataPrefix = self.config.coaddName + "Coadd_"

    @classmethod
    def _makeArgumentParser(cls):
        parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "deepCoadd", help="data ID, e.g. --id tract=12345 patch=1,2",
                               ContainerClass=CoaddDataIdContainer)
        return parser

    def getPreviousTaskClass(self):
        return MeasureCoaddTask

    def _getConfigName(self):
        """Return the name of the config dataset
        """
        return "%s_measureCoadd_config" % (self.config.coaddName,)

    def _getMetadataName(self):
        """Return the name of the metadata dataset
        """
        return "%s_measureCoadd_metadata" % (self.config.coaddName,)
