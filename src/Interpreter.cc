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

#include "lsst/meas/multifit/Interpreter.h"

namespace lsst { namespace meas { namespace multifit {

Interpreter::Interpreter(int parameterDim, PTR(Model) model, PTR(Prior) prior)
  : _parameterDim(parameterDim), _model(model), _prior(prior)
{}

void Interpreter::packParameters(
    ndarray::Array<Scalar const,1,1> const & nonlinear,
    ndarray::Array<Scalar const,1,1> const & amplitudes,
    ndarray::Array<Scalar,1,1> const & parameters
) const {
    LSST_ASSERT_EQUAL(
        parameters.getSize<0>(), getParameterDim(),
        "Size of parameter array (%d) does not match expected size (%d)",
        pex::exceptions::LengthErrorException
    );
    LSST_ASSERT_EQUAL(
        nonlinear.getSize<0>(), getNonlinearDim(),
        "Size of nonlinear array (%d) does not match expected size (%d)",
        pex::exceptions::LengthErrorException
    );
    LSST_ASSERT_EQUAL(
        amplitudes.getSize<0>(), getAmplitudeDim(),
        "Size of ampiltude array (%d) does not match expected size (%d)",
        pex::exceptions::LengthErrorException
    );
    _packParameters(nonlinear, amplitudes, parameters);
}

void Interpreter::unpackNonlinear(
    ndarray::Array<Scalar const,1,1> const & parameters,
    ndarray::Array<Scalar,1,1> const & nonlinear
) const {
    LSST_ASSERT_EQUAL(
        parameters.getSize<0>(), getParameterDim(),
        "Size of parameter array (%d) does not match expected size (%d)",
        pex::exceptions::LengthErrorException
    );
    LSST_ASSERT_EQUAL(
        nonlinear.getSize<0>(), getNonlinearDim(),
        "Size of nonlinear array (%d) does not match expected size (%d)",
        pex::exceptions::LengthErrorException
    );
    _unpackNonlinear(parameters, nonlinear);
}

ndarray::Array<Scalar,2,2> Interpreter::computeParameterQuantiles(
    ModelFitRecord const & record,
    ndarray::Array<Scalar const,1,1> const & fractions
) const {
    ndarray::Array<Scalar,2,2> output = ndarray::allocate(getParameterDim(), fractions.getSize<0>());
    for (int i = 0; i < getParameterDim(); ++i) {
        output[i] = computeParameterQuantiles(record, fractions, i);
    }
    return output;
}

ndarray::Array<Scalar,2,2> Interpreter::computeNonlinearQuantiles(
    ModelFitRecord const & record,
    ndarray::Array<Scalar const,1,1> const & fractions
) const {
    ndarray::Array<Scalar,2,2> output = ndarray::allocate(getNonlinearDim(), fractions.getSize<0>());
    for (int i = 0; i < getNonlinearDim(); ++i) {
        output[i] = computeNonlinearQuantiles(record, fractions, i);
    }
    return output;
}

ndarray::Array<Scalar,2,2> Interpreter::computeAmplitudeQuantiles(
    ModelFitRecord const & record,
    ndarray::Array<Scalar const,1,1> const & fractions
) const {
    ndarray::Array<Scalar,2,2> output = ndarray::allocate(getAmplitudeDim(), fractions.getSize<0>());
    for (int i = 0; i < getAmplitudeDim(); ++i) {
        output[i] = computeAmplitudeQuantiles(record, fractions, i);
    }
    return output;
}

ndarray::Array<Scalar,2,2> Interpreter::computeParameterCovariance(ModelFitRecord const & record) const {
    return computeParameterCovariance(record, computeParameterMean(record));
}

ndarray::Array<Scalar,2,2> Interpreter::computeNonlinearCovariance(ModelFitRecord const & record) const {
    return computeNonlinearCovariance(record, computeNonlinearMean(record));
}

ndarray::Array<Scalar,2,2> Interpreter::computeAmplitudeCovariance(ModelFitRecord const & record) const {
    return computeAmplitudeCovariance(record, computeAmplitudeMean(record));
}


}}} // namespace lsst::meas::multifit
