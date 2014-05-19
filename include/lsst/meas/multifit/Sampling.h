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

#ifndef LSST_MEAS_MULTIFIT_Sampling_h_INCLUDED
#define LSST_MEAS_MULTIFIT_Sampling_h_INCLUDED

#include "ndarray.h"
#include "lsst/afw/table/fwd.h"
#include "lsst/afw/table/Key.h"
#include "lsst/meas/multifit/Mixture.h"
#include "lsst/meas/multifit/Likelihood.h"
#include "lsst/meas/multifit/Interpreter.h"

namespace lsst { namespace meas { namespace multifit {

class SamplingObjective;

/**
 *  @brief Interpreter for Monte Carlo sampling
 *
 *  In addition to serving as objects that can be used to interpret fitting result after the fit,
 *  SamplingInterpreters also serve as factories for SamplingObjectives (via a friend function,
 *  makeSamplingObjective).
 *
 *  This class is an intermediate base class for the two concrete interpreters (DirectSamplingInterpreter
 *  and MarginalSamplingInterpreter) that represent different ways of handling the amplitude parameters
 *  of the model.
 */
class SamplingInterpreter : public Interpreter {
public:

    ArrayKey const & getParameterKey() const { return _parameterKey; }

    ArrayKey const & getNonlinearKey() const { return _nonlinearKey; }

    ScalarKey const & getWeightKey() const { return _weightKey; }

    virtual ndarray::Array<Scalar,1,1> computeParameterQuantiles(
        ModelFitRecord const & record,
        ndarray::Array<Scalar const,1,1> const & fractions,
        int index
    ) const;

    virtual ndarray::Array<Scalar,1,1> computeNonlinearQuantiles(
        ModelFitRecord const & record,
        ndarray::Array<Scalar const,1,1> const & fractions,
        int index
    ) const;

    virtual ndarray::Array<Scalar,1,1> computeParameterMean(
        ModelFitRecord const & record
    ) const;

    virtual ndarray::Array<Scalar,1,1> computeNonlinearMean(
        ModelFitRecord const & record
    ) const;

    virtual ndarray::Array<Scalar,2,2> computeParameterCovariance(
        ModelFitRecord const & record,
        ndarray::Array<Scalar const,1,1> const & mean
    ) const;

    virtual ndarray::Array<Scalar,2,2> computeNonlinearCovariance(
        ModelFitRecord const & record,
        ndarray::Array<Scalar const,1,1> const & mean
    ) const;

    friend PTR(SamplingObjective) makeSamplingObjective(
        PTR(SamplingInterpreter) interpreter,
        PTR(Likelihood) likelihood
    );

protected:

    // Protected constructor; subclasses must also initialize _parameterKey and _nestedKey in their ctor body
    // (it just works better that way)
    SamplingInterpreter(
        afw::table::Schema & sampleSchema,
        Model::NameVector const & parameterNames,
        PTR(Model) model,
        PTR(Prior) prior
    );

    virtual PTR(SamplingObjective) makeObjective(
        PTR(SamplingInterpreter) self,
        PTR(Likelihood) likelihood
    ) const = 0;

    ndarray::Array<Scalar,1,1> computeSampleQuantiles(
        afw::table::BaseCatalog const & samples,
        ndarray::Array<Scalar const,1,1> const & fractions,
        ScalarKey const & targetKey
    ) const;

    ndarray::Array<Scalar,1,1> computeSampleMean(
        afw::table::BaseCatalog const & samples,
        ArrayKey const & targetKey
    ) const;

    ndarray::Array<Scalar,2,2> computeSampleCovariance(
        afw::table::BaseCatalog const & samples,
        ndarray::Array<Scalar const,1,1> const & mean,
        ArrayKey const & targetKey
    ) const;

    ArrayKey _parameterKey;
    ArrayKey _nonlinearKey;
    ScalarKey _weightKey;
};

class SamplingObjective {
public:

    int getParameterDim() const { return _interpreter->getParameterDim(); }

    PTR(SamplingInterpreter) getInterpreter() const { return _interpreter; }

    virtual Scalar operator()(
        ndarray::Array<Scalar const,1,1> const & parameters,
        afw::table::BaseRecord & sample
    ) const = 0;

    virtual ~SamplingObjective() {}

protected:
    SamplingObjective(PTR(SamplingInterpreter) interpreter, PTR(Likelihood) likelihood);

    PTR(SamplingInterpreter) _interpreter;
    PTR(Likelihood) _likelihood;
    ndarray::Array<Pixel,2,-1> _modelMatrix;
};

/**
 *  @brief Factory function used to create SamplingObjectives
 *
 *  If the given Prior is not null, it will be used instead of any Prior attached to the
 *  Interpreter.
 *
 *  The actual work is delegated to a protected member function of the SamplingInterpreter,
 *  but we have to provide this interface so we can provide it with a shared_ptr to itself
 *  (this seems preferable to making the entire hierarchy inherit from enable_shared_from_this).
 */
inline PTR(SamplingObjective) makeSamplingObjective(
    PTR(SamplingInterpreter) interpreter,
    PTR(Likelihood) likelihood
) {
    return interpreter->makeObjective(interpreter, likelihood);
}

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_Sampling_h_INCLUDED
