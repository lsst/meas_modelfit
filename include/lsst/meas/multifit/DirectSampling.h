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

#ifndef LSST_MEAS_MULTIFIT_DirectSampling_h_INCLUDED
#define LSST_MEAS_MULTIFIT_DirectSampling_h_INCLUDED

#include "ndarray.h"
#include "lsst/meas/multifit/Sampling.h"

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief A SamplingInterpreter when we don't marginalize amplitudes
 *
 *  In this case, the parameter vector is just the concatenation of the nonlinear and amplitude vectors
 *  (in that order).  It's all quite straightforward, and its easy to interpret the results - the only
 *  downside is that we have to explore a higher-dimensional parameter space.
 */
class DirectSamplingInterpreter : public SamplingInterpreter {
public:

    DirectSamplingInterpreter(
        afw::table::Schema & sampleSchema,
        PTR(Model) model,
        PTR(Prior) prior=PTR(Prior)()
    );

    PTR(DirectSamplingInterpreter) clone() const {
        return boost::static_pointer_cast<DirectSamplingInterpreter>(_clone());
    }

    virtual ndarray::Array<Scalar,1,1> computeAmplitudeQuantiles(
        ModelFitRecord const & record,
        ndarray::Array<Scalar const,1,1> const & fractions,
        int index
    ) const;

    virtual ndarray::Array<Scalar,1,1> computeAmplitudeMean(
        ModelFitRecord const & record
    ) const;

    virtual ndarray::Array<Scalar,2,2> computeAmplitudeCovariance(
        ModelFitRecord const & record,
        ndarray::Array<Scalar const,1,1> const & mean
    ) const;

protected:

    virtual PTR(Interpreter) _clone() const;

    virtual PTR(SamplingObjective) makeObjective(
        PTR(SamplingInterpreter) self,
        PTR(Likelihood) likelihood
    ) const;

    virtual void _packParameters(
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        ndarray::Array<Scalar const,1,1> const & amplitudes,
        ndarray::Array<Scalar,1,1> const & parameters
    ) const;

    virtual void _unpackNonlinear(
        ndarray::Array<Scalar const,1,1> const & parameters,
        ndarray::Array<Scalar,1,1> const & nonlinear
    ) const;

};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_MarginalSampling_h_INCLUDED
