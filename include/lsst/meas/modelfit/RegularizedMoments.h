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

#include <cmath>
#include <exception>
#include <string>
#include <tuple>
#include <vector>
#include <Eigen/Dense>

namespace lsst {
namespace meas {
namespace modelfit {

struct Moments {
    using FirstMoment = Eigen::Matrix<double, 2, 1>
    using SecondMoment = Eigen::Matrix<double, 2, 2>
    using ParameterVector = Eigen::Matrix<double, 6, 1>
    using ParameterMatrix = Eigen::Matrix<double, 6, 6>

    template <typename iter>
    Moments(iter begin, iter end);

    Moments(std::vector<double> moments);

    Moments(std::initializer_list<double> initList);

    Moments(double zero, FirstMoment first, SecondMoment second);

    double operator[](int i);

    ParameterVector getParameterVector();

    bool aproxEqual(Moments const & other, double tol=1e-6);

    double zeroth;
    FirstMoment first;
    SecondMoment second;
}

struct ZerothShapeMoment {
    static double computeValue(Moments Q, Moments W);

    static Moments::ParameterVector computeGradient(Moments Q, Moments W);
};

struct FirstShapeMoment {
    static double computeValue(Moments Q, Moments W);

    static Moments::ParameterVector computeGradient(Moments Q, Moments W);
};

struct SecondShapeMoment {
    static double computeValue(Moments Q, Moments W);

    static Moments::ParameterVector computeGradient(Moments Q, Moments W);
};

// Tests for classes in anonymous name spaces
std::pair<Moments, Moments> buildTestMoments();
bool testAlphaX(double tol=1e-6);
}}} // Close namespace lsst::meas::modelfit
