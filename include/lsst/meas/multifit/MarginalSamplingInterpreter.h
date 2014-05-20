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

#ifndef LSST_MEAS_MULTIFIT_MarginalSamplingInterpreter_h_INCLUDED
#define LSST_MEAS_MULTIFIT_MarginalSamplingInterpreter_h_INCLUDED

#include "ndarray.h"
#include "lsst/afw/math/Random.h"
#include "lsst/meas/multifit/Sampling.h"
#include "lsst/meas/multifit/ModelFitRecord.h"

namespace lsst { namespace meas { namespace multifit {

class DirectSamplingInterpreter;

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
 *  Right now, the nested information isn't used, and we just return NaNs when someone tries to get
 *  amplitude information out of a marginalized distribution.  Note that the difficulty isn't unpacking
 *  the (Gaussian) likelihood, but rather combining that with the (non-Gaussian) prior.
 */
class MarginalSamplingInterpreter : public SamplingInterpreter {
public:

    MarginalSamplingInterpreter(
        afw::table::Schema & sampleSchema,
        PTR(Model) model,
        PTR(Prior) prior=PTR(Prior)()
    );

    ArrayKey getNestedKey() const { return _nestedKey; }

#ifndef SWIG // can't use Eigen references as Python output arguments

    void unpackNested(
        ndarray::Array<Scalar const,1,1> const & nested, Vector & gradient, Matrix & hessian
    ) const;

    void unpackNested(
        afw::table::BaseRecord const & sample, Vector & gradient, Matrix & hessian
    ) const {
        unpackNested(sample.get(_nestedKey), gradient, hessian);
    }

#endif

    void unpackNested(
        ndarray::Array<Scalar const,1,1> const & nested,
        ndarray::Array<Scalar,1,1> const & gradient,
        ndarray::Array<Scalar,2,2> const & hessian
    ) const;

    void unpackNested(
        afw::table::BaseRecord const & sample,
        ndarray::Array<Scalar,1,1> const & gradient,
        ndarray::Array<Scalar,2,2> const & hessian
    ) const {
        unpackNested(sample.get(_nestedKey), gradient, hessian);
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

/**
 *  @brief Convert the samples and proposal pdf done with marginal sampling to direct sampling.
 *
 *  This functor is initialized with the marginal and direct Schemas and Interpreters, and uses
 *  these to draw new samples that include the amplitudes in the parameters, and then create a
 *  new proposal that includes amplitude dimensions that match those samples.
 *
 *  This object does not handle mapping fields in the ModelFitRecords themselves (that's a task
 *  left for Python code, which handles those Schemas).
 */
class UnnestMarginalSamples {
public:

    /**
     *  @brief Construct from existing schemas and interpreters
     *
     *  If the multiplier is <= 0, it will be set to the ratio of the number of direct parameters
     *  to marginal parameters (i.e. (nonlinearDim + amplitudeDim) / nonlinearDim).
     */
    UnnestMarginalSamples(
        afw::table::Schema const & marginalSchema,
        afw::table::Schema const & directSchema,
        PTR(MarginalSamplingInterpreter) marginalInterpreter,
        PTR(DirectSamplingInterpreter) directInterpreter,
        PTR(afw::math::Random) rng,
        double multiplier=-1.0
    );

    /// Create a new sample catalog and proposal pdf from marginalRecord, and add them to directRecord.
    void apply(ModelFitRecord const & marginalRecord, ModelFitRecord & directRecord) const;

private:
    afw::table::SchemaMapper _mapper;
    double _multiplier;
    PTR(afw::math::Random) _rng;
    PTR(MarginalSamplingInterpreter) _marginalInterpreter;
    PTR(DirectSamplingInterpreter) _directInterpreter;
    PTR(Prior) _prior;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_MarginalSamplingInterpreter_h_INCLUDED
