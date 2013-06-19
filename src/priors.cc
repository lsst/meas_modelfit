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

#include "boost/math/special_functions/erf.hpp"
#include "Eigen/LU"

#include "lsst/pex/exceptions.h"
#include "lsst/meas/multifit/priors.h"

namespace lsst { namespace meas { namespace multifit {

namespace {

double solveFlat(LogGaussian const & likelihood, samples::Vector & z) {
    static double const LOG_2PI = std::log(2.0 * M_PI);
    Eigen::FullPivLU<samples::Matrix> lu(likelihood.fisher);
    if (!lu.isInvertible()) {
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            "Amplitude Fisher matrix is singular; please use a non-flat Prior"
        );
    }
    z = lu.solve(likelihood.grad);
    return - 0.5 * likelihood.grad.dot(z) - 0.5*LOG_2PI*likelihood.grad.size()
        + 0.5*std::log(lu.determinant());
}

} // anonymous

PTR(FlatPrior) FlatPrior::get() {
    static PTR(FlatPrior) instance(new FlatPrior());
    return instance;
}

samples::Scalar FlatPrior::apply(LogGaussian const & likelihood, samples::Vector const & parameters) const {
    samples::Vector z(likelihood.grad.size());
    return solveFlat(likelihood, z);
}

samples::Scalar FlatPrior::computeFluxExpectation(
    LogGaussian const & likelihood, samples::Vector const & parameters
) const {
    samples::Vector z(likelihood.grad.size());
    solveFlat(likelihood, z);
    return -z.sum();
}

samples::Scalar FlatPrior::computeSquaredFluxExpectation(
    LogGaussian const & likelihood, samples::Vector const & parameters
) const {
    throw LSST_EXCEPT(
        pex::exceptions::LogicErrorException,
        "NOT IMPLEMENTED"
    );
}

samples::Vector FlatPrior::computeFractionExpectation(
    LogGaussian const & likelihood, samples::Vector const & parameters
) const {
    throw LSST_EXCEPT(
        pex::exceptions::LogicErrorException,
        "NOT IMPLEMENTED"
    );
}

}}} // namespace lsst::meas::multifit
