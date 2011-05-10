// -*- lsst-c++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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

%{
#include "lsst/meas/multifit/BaseInterpreter.h"
#include "lsst/meas/multifit/SimpleInterpreter.h"
%}

%declareNumPyConverters(Eigen::Matrix5d);
%declareNumPyConverters(Eigen::Matrix2d);

SWIG_SHARED_PTR(BaseInterpreterPtr, lsst::meas::multifit::BaseInterpreter);
SWIG_SHARED_PTR_DERIVED(UnifiedInterpreterPtr, 
        lsst::meas::multifit::BaseInterpreter,
        lsst::meas::multifit::UnifiedInterpreter);
SWIG_SHARED_PTR_DERIVED(NestedInterpreterPtr, 
        lsst::meas::multifit::BaseInterpreter,
        lsst::meas::multifit::NestedInterpreter);

SWIG_SHARED_PTR_DERIVED(SimpleInterpreterPtr
        lsst::meas::multifit::BaseInterpreter
        lsst::meas::multifit::SimpleInterpreter);
SWIG_SHARED_PTR_DERIVED(UnifiedSimpleInterpreterPtr, 
        lsst::meas::multifit::SimpleInterpreter,
        lsst::meas::multifit::UnifiedInterpreter);
SWIG_SHARED_PTR_DERIVED(NestedSimpleInterpreterPtr, 
        lsst::meas::multifit::SimpleInterpreter,
        lsst::meas::multifit::NestedSimpleInterpreter);


%include "lsst/meas/multifit/BaseInterpreter.h"
%include "lsst/meas/multifit/SimpleInterpreter.h"
