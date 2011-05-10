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
#include "lsst/meas/multifit/BaseDistribution.h"
#include "lsst/meas/multifit/SimpleDistribution.h"
#include "lsst/meas/multifit/StudentDistribution.h"
#include "lsst/meas/multifit/GaussianDistribution.h"
%}

%declareNumPyConverters(lsst::ndarray::Array<double, 1, 1>);
%declareNumPyConverters(lsst::ndarray::Array<double const, 1, 1>);
%declareNumPyConverters(lsst::ndarray::Array<double const, 2, 1>);

SWIG_SHARED_PTR(BaseDistributionPtr, lsst::meas::multifit::BaseDistribution);
SWIG_SHARED_PTR_DERIVED(SimpleDistributionPtr, 
        lsst::meas::multifit::BaseDistribution, 
        lsst::meas::multifit::SimpleDistribution);
SWIG_SHARED_PTR_DERIVED(StudentDistributionPtr, 
        lsst::meas::multifit::SimpleDistribution
        lsst::meas::multifit::StudentDistribution);
SWIG_SHARED_PTR_DERIVED(GaussianDistributionPtr, 
        lsst::meas::multifit::SimpleDistribution
        lsst::meas::multifit::GaussianDistribution);

%include "lsst/meas/multifit/BaseDistribution.h"
%include "lsst/meas/multifit/SimpleDistribution.h"
%include "lsst/meas/multifit/StudentDistribution.h"
%include "lsst/meas/multifit/GaussianDistribution.h"
