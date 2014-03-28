// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

%shared_ptr(lsst::meas::multifit::CModelControl)
%shared_ptr(lsst::meas::multifit::CModelAlgorithm)
%include "lsst/meas/multifit/CModel.h"

%extend lsst::meas::multifit::CModelStageResult {
%pythoncode %{

def displayHistory(self, *kwds):
    """Return a display.xOptimizerDisplay object that shows the track of the optimizer
    in this fit.  Additional keyword arguments are forwarded to the OptimizerDisplay
    constructor.
    """
    from .display import OptimizerDisplay
    return OptimizerDisplay(self.history, self.model, self.objfunc, *kwds)

%}
}

%pythoncode %{
import lsst.pex.config
import lsst.meas.algorithms

CModelStageConfig = lsst.pex.config.makeConfigClass(CModelStageControl)

CModelRegionConfig = lsst.pex.config.makeConfigClass(CModelRegionControl)

CModelDiagnosticsConfig = lsst.pex.config.makeConfigClass(CModelDiagnosticsControl)

lsst.meas.algorithms.AlgorithmRegistry.register("cmodel", CModelControl)

CModelAlgorithm.Result = CModelResult

%}
