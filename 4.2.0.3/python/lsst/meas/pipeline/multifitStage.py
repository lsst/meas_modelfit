#! /usr/bin/env python

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


import lsst.pex.harness.stage as harnessStage
from lsst.pex.logging import Log
import lsst.daf.base as dafBase
from lsst.daf.base import *
import lsst.pex.policy as pexPolicy
import lsst.meas.multifit as measMult
import lsst.meas.utils.multifit as utilsMult

__all__ = ["MultifitStage", "MultifitStageParallel"]

class MultifitStageParallel(harnessStage.ParallelProcessing):
    """
    Given an Exposure Stack and an initial Model, fit the model on the
    stack using multifit.

    INPUT:
    - initial model
    - ExposureList
    OUTPUT:
    - fit model
    - double : sgChisq
    - double : psChisq
    """
    
    def setup(self):
        self.log = Log(self.log, "MultifitStage - parallel")
        policyFile = pexPolicy.DefaultPolicyFile("meas_pipeline", 
            "MultifitStageDictionary.paf", "policy")
        defPolicy = pexPolicy.Policy.createPolicy(policyFile, policyFile.getRepositoryPah(), True)

        if self.policy is None:
            self.policy = pexPolicy.Policy()
        self.policy.mergeDefaults(defPolicy.getDictionary())

    def process(self, clipboard):
        exposureList = clipboard.get(self.policy.get("input.exposureList"))

        fitter = measMult.SingleLinearParameterFitter(
            self.policy.get("parameters.fitterPolicy")
        )

        model = clipboard.get(self.policy.get("input.model"))
        evaluator = measMult.ModelEvaluator(model)
        evaluator.setExposureList(exposureList)
        psResult = fitter.apply(evaluator)
        clipboard.put(self.policy.get("output.model"), psResult.getModel())
        clipboard.put(self.policy.get("output.chisq"), psResult.getChisq())

class MultifitStage(harnessStage.Stage):
    parallelClass = MultifitStageParallel
