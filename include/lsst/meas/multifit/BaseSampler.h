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

#ifndef LSST_MEAS_MULTIFIT_BaseSampler_h_INCLUDED
#define LSST_MEAS_MULTIFIT_BaseSampler_h_INCLUDED

#include "lsst/meas/multifit/SampleSet.h"
#include "lsst/meas/multifit/Likelihood.h"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief C++ base class for likelihood samplers
 *
 *  The "Sampler" classes are actually tasks defined in Python, but the real (non-bookkeeping) work
 *  is generally done here.  There should be one BaseSampler subclass for each BaseSamplerTask
 *  subclass, and its up to them to determine how to construct the state object in the Task's setup()
 *  and reset() methods.  BaseSamplerTask.run() will then delegate to BaseSamplerState::run() to do
 *  the actual work of filling a SampleSet with new likelihood samples.
 *
 *  Note that there is one BaseSampler instance per object, but one SamplerTask can be run on many
 *  objects, and hence its attributes do not include any per-object state.  That's the job of this class,
 *  which is created by some of the Task's methods and passed to others.
 */
class BaseSampler {
public:

    /// Draw and evaluate samples using the given object.
    virtual SampleSet run(Likelihood const & likelihood) const = 0;

    virtual ~BaseSampler() {}

};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_BaseSampler_h_INCLUDED
