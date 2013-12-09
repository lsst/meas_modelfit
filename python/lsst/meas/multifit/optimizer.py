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

import numpy

import lsst.pex.config
import lsst.pipe.base
import lsst.afw.table

from . import multifitLib

__all__ = ("OptimizerConfig", "OptimizerTask")

@lsst.pex.config.wrap(multifitLib.OptimizerControl)
class OptimizerConfig(lsst.pex.config.Config):
    doRecordHistory = lsst.pex.config.Field(
        dtype=bool, default=False,
        doc="Whether to save optimizer tracks in the samples table"
        )
    doRecordDerivatives = lsst.pex.config.Field(
        dtype=bool, default=True,
        doc="Whether to save derivatives with history (ignored if doRecordHistory is False)"
        )

class OptimizerTask(lsst.pipe.base.Task):
    """A 'fitter' subtask for Measure tasks that uses a greedy optimizer.
    """

    ConfigClass = OptimizerConfig

    def __init__(self, schema, keys, model, prior, **kwds):
        lsst.pipe.base.Task.__init__(self, **kwds)
        # n.b. schema argument is for modelfits catalog; self.sampleSchema is for sample catalog
        self.sampleSchema = lsst.afw.table.Schema()
        self.interpreter = multifitLib.OptimizerInterpreter(model, prior)
        self.keys = keys
        self.keys["ref.parameters"] = schema.addField(
            "ref.parameters", type="ArrayD", size=self.interpreter.getParameterDim(),
            doc="sampler parameters from reference catalog"
            )
        self.keys["fit.parameters"] = schema.addField(
            "fit.parameters", type="ArrayD", size=self.interpreter.getParameterDim(),
            doc="best-fit sampler parameters"
            )
        self.keys["fit.flags"] = schema.addField("fit.flags", type="Flag", doc="whether the fit converged")
        # TODO: should use true Flag fields in production, but that's a lot of boilerplate - we should have
        # a class for bitflags that we can convert to table Flags in bulk
        self.keys["fit.state"] = schema.addField("fit.state", type=int,
                                                 doc="State flags transferred directly from Optimizer")
        self.keys["fit.objective"] = schema.addField("objective", type=float,
                                                     doc="objective (-ln likelihood*prior) at fit.parameters")
        if self.config.doRecordHistory:
            self.recorder = multifitLib.OptimizerHistoryRecorder(
                self.sampleSchema, model, self.config.doRecordDerivatives
                )
        else:
            self.recorder = None

    def makeTable(self):
        """Return a Table object that can be used to construct sample records.
        """
        return lsst.afw.table.BaseTable.make(self.sampleSchema)

    def initialize(self, record):
        """Initialize an output record, setting any derived fields and record
        attributes (i.e. samples or pdf) needed before calling run().

        This method is not called when using a "warm start" from a previous fit.
        """
        parameters = record[self.keys["ref.parameters"]]
        self.interpreter.packParameters(record[self.keys["ref.nonlinear"]],
                                        record[self.keys["ref.amplitudes"]],
                                        parameters)
        record[self.keys["fit.parameters"]] = parameters

    def run(self, likelihood, record):
        """Do the actual fitting, using the given likelihood, update the 'pdf' and 'samples' attributes,
        and save best-fit values in the 'fit.parameters' field.
        """
        parameters = record[self.keys["fit.parameters"]]
        objective = multifitLib.OptimizerObjective.makeFromLikelihood(likelihood, self.interpreter.getPrior())
        optimizer = multifitLib.Optimizer(objective, parameters, self.config.makeControl())
        if self.recorder:
            optimizer.run(self.recorder, record.getSamples())
        else:
            optimizer.run()
        self.interpreter.attachPdf(record, optimizer)
        parameters[:] = self.interpreter.computeParameterMean(record)
        record.set(self.keys['fit.flags'], bool(optimizer.getState() & multifitLib.Optimizer.FAILED))
        record.set(self.keys['fit.state'], optimizer.getState())
        record.set(self.keys['fit.objective'], optimizer.getObjectiveValue())
