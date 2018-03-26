// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2018  AURA/LSST.
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
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */
#ifndef LSST_MEAS_MODELFIT_REGULARIZEDMOMENTS_H
#define LSST_MEAS_MODELFIT_REGULARIZEDMOMENTS_H

#include <Eigen/Dense>

namespace lsst {
namespace meas {
namespace modelfit {

class MomentsModel {
public:
    using Moments = Eigen::Matrix<double, 6, 1>;
    using Jacobian = Eigen::Matrix<double, 6, 6>;
    using FirstMoment = Eigen::Matrix<double, 2, 1>;
    using SecondMoment = Eigen::Matrix<double, 2, 2>;

    MomentsModel(Moments const & W): W(W) {
    }

    void at(Moments const & Q);

    Moments computeValues();

    Jacobian computeJacobian();

private:
    void makeValue();

    Moments W, Q, value;

    double norm;

    FirstMoment alpha;

    SecondMoment beta;
};

// Tests for classes in anonymous name spaces
bool testNorm(double tol=1e-6);

bool testAlphaX(double tol=1e-6);
bool testAlphaY(double tol=1e-6);

bool testBetaX(double tol=1e-6);
bool testBetaY(double tol=1e-6);
bool testBetaXY(double tol=1e-6);
}}} // Close namespace lsst::meas::modelfit

#endif // LSST_MEAS_MODELFIT_REGULARIZEDMOMENTS_H

