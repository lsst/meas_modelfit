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

#ifndef LSST_MEAS_MULTIFIT_MarginalSampling_h_INCLUDED
#define LSST_MEAS_MULTIFIT_MarginalSampling_h_INCLUDED

#include "ndarray.h"
#include "lsst/meas/multifit/Sampling.h"

namespace lsst { namespace meas { namespace multifit {


/**
 *  @brief A SamplingInterpreter for when we marginalize over amplitudes at each parameter point.
 *
 *  In this case, the parameter vector is just the nonlinear vector, and we store information about
 *  the amplitude likelihood in another "nested" field we add to the sample schema.  That "nested" data
 *  is intended to be opaque, but it's just the gradient and Hessian of the amplitude likelihood packed
 *  into a 1-d array.  Using NumPy syntax, that's (for an N-d amplitude):
 *  @code
 *    nested[:N] = gradient
 *    nested[N] = hessian[0,0]
 *    nested[N:N+1] = hessian[1,:1]
 *    nested[N+1:N+3] = hessian[2,:2]
 *    ...
 *  @endcode
 *  That is, after the gradient, the symmetric Hessian is stored packed in what LAPACK calls UPLO='U' order.
 *
 *  Right now, the nested information isn't used, and we just throw exceptions when someone tries to get
 *  amplitude information out of a marginalized distribution.  Note that the difficulty isn't unpacking
 *  the (Gaussian) likelihood, but rather combining that with the (non-Gaussian) prior.
 */
class MarginalSamplingInterpreter : public SamplingInterpreter {
public:

    MarginalSamplingInterpreter(
        afw::table::Schema & sampleSchema,
        PTR(Model) model,
        PTR(Prior) prior
    );

    ArrayKey getNestedKey() const { return _nestedKey; }

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

private:
    ArrayKey _nestedKey;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_MarginalSampling_h_INCLUDED
