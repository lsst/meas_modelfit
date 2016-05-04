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

%shared_ptr(lsst::meas::modelfit::CModelControl)
%shared_ptr(lsst::meas::modelfit::CModelAlgorithm)
%newobject lsst::meas::modelfit::CModelAlgorithm::apply;
%newobject lsst::meas::modelfit::CModelAlgorithm::applyForced;

// SWIG accessing this directly is a bad idea: it's a polymorphic type that would be returned by value,
// but SWIG wants to do everything by pointer.  The result is confusion over lifetimes, and dangling pointers.
%ignore lsst::meas::modelfit::CModelStageResult::ellipse;

%include "lsst/meas/modelfit/PixelFitRegion.h"
%include "lsst/meas/modelfit/CModel.h"

%extend lsst::meas::modelfit::CModelStageResult {

    // An alternative accessor for ellipse
    std::shared_ptr<lsst::afw::geom::ellipses::Quadrupole> getEllipse() const {
        return std::make_shared<lsst::afw::geom::ellipses::Quadrupole>($self->ellipse);
    }

%pythoncode %{

def displayHistory(self, *kwds):
    """Return a display.OptimizerDisplay object that shows the track of the optimizer
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

PixelFitRegionConfig = lsst.pex.config.makeConfigClass(PixelFitRegionControl)

CModelConfig = lsst.pex.config.makeConfigClass(CModelControl)

CModelAlgorithm.Result = CModelResult
CModelAlgorithm.ConfigClass = CModelConfig

%}
