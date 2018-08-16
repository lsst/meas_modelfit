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
"""
A script to be used when profiling and debugging large-catalog IO.
"""
import resource
import os
import lsst.meas.modelfit


def main(nRecords, nSamplesPerRecord):
    print("Initializing Task")
    task = lsst.meas.modelfit.MeasureCcdTask()
    print("Filling catalog")
    catalog = lsst.meas.modelfit.ModelFitCatalog(task.makeTable())
    for i in range(nRecords):
        record = catalog.addNew()
        samples = record.getSamples()
        for j in range(nSamplesPerRecord):
            samples.addNew()
    print("Saving catalog")
    res0 = resource.getrusage(resource.RUSAGE_SELF)
    catalog.writeFits("tmp.fits")
    res1 = resource.getrusage(resource.RUSAGE_SELF)
    print("Save complete: system=%f, user=%f" %
          (res1.ru_stime - res0.ru_utime, res1.ru_utime - res0.ru_utime))
    print("Loading catalog")
    res0 = resource.getrusage(resource.RUSAGE_SELF)
    lsst.meas.modelfit.ModelFitCatalog.readFits("tmp.fits")
    res1 = resource.getrusage(resource.RUSAGE_SELF)
    print("Load complete: system=%f, user=%f" %
          (res1.ru_stime - res0.ru_utime, res1.ru_utime - res0.ru_utime))
    os.remove("tmp.fits")


if __name__ == "__main__":
    import sys
    main(int(sys.argv[1]), int(sys.argv[2]))
