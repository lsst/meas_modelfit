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

#ifndef LSST_MEAS_MULTIFIT_Sampler_h_INCLUDED
#define LSST_MEAS_MULTIFIT_Sampler_h_INCLUDED

#include "ndarray.h"
#include "lsst/afw/table/fwd.h"
#include "lsst/afw/table/Key.h"
#include "lsst/meas/multifit/models.h"
#include "lsst/meas/multifit/Mixture.h"
#include "lsst/meas/multifit/Likelihood.h"
#include "lsst/meas/multifit/priors.h"

namespace lsst { namespace meas { namespace multifit {

class SamplerObjective {
public:

    typedef afw::table::Key< afw::table::Array<Scalar> > ArrayKey;

    int getParameterDim() const { return _parameterKey.getSize(); }

    ArrayKey const & getParameterKey() const { return _parameterKey; }

    virtual Scalar operator()(
        ndarray::Array<Scalar const,1,1> const & parameters,
        afw::table::BaseRecord & sample
    ) const = 0;

    virtual ~SamplerObjective() {}

protected:
    SamplerObjective(ArrayKey const & parameterKey, PTR(Likelihood) likelihood, PTR(Prior) prior);

    ArrayKey _parameterKey;
    PTR(Likelihood) _likelihood;
    PTR(Prior) _prior;
    ndarray::Array<Pixel,2,-1> _modelMatrix;
};

class SamplerObjectiveFactory {
public:

    PTR(SamplerObjective) operator()(PTR(Likelihood) likelihood) const;

    int getParameterDim() const { return _parameterKey.getSize(); }

    SamplerObjective::ArrayKey const & getParameterKey() const { return _parameterKey; }

    void mapParameters(
        ndarray::Array<Scalar const,1,1> const & nonlinear,
        ndarray::Array<Scalar const,1,1> const & amplitudes,
        ndarray::Array<Scalar,1,1> const & parameters
    ) const;

    SamplerObjectiveFactory(
        afw::table::Schema & sampleSchema,
        PTR(Model) model,
        PTR(Prior) prior,
        bool doMarginalizeAmplitudes=true
    );

private:
    bool _doMarginalizeAmplitudes;
    PTR(Prior) _prior;
    SamplerObjective::ArrayKey _parameterKey;
    SamplerObjective::ArrayKey _nestedKey;
};

class Sampler {
public:

    virtual void run(
        SamplerObjective const & objective,
        PTR(Mixture) proposal,
        afw::table::BaseCatalog & samples
    ) const = 0;

    virtual ~Sampler() {}

};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_Sampler_h_INCLUDED
