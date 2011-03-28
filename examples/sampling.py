#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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

import lsst.afw.image
import lsst.afw.detection
import lsst.afw.image
import lsst.afw.geom.ellipses
import lsst.meas.multifit.viewer
import lsst.meas.multifit.sampling
import numpy
from matplotlib import pyplot

def main(viewer, size=10000):
    engine = lsst.meas.multifit.sampling.RandomEngine()
    mean = viewer.parameters.copy()
    sigma = numpy.identity(mean.size, dtype=float)
    sigma[0,0] = 3.0
    sigma[1,1] = 0.3
    sigma[2,2] = 0.3
    importance = lsst.meas.multifit.sampling.MixtureDistribution(
        [lsst.meas.multifit.sampling.MixtureComponent(1.0, mean, sigma)], -1,
        )
    sampler = lsst.meas.multifit.sampling.IterativeImportanceSampler(
        viewer.evaluator, importance, engine
        )
    sampler.run(size)
    return sampler
