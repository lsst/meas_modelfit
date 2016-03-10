// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2015 LSST/AURA
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

#include "Eigen/LU"

#include "lsst/meas/modelfit/detail/polynomials.h"

namespace lsst { namespace meas { namespace modelfit { namespace detail {


template <int N>
typename Vandermonde<N>::RowVector
Vandermonde<N>::eval(double x) {
    RowVector z = RowVector::Zero();
    double y = 1.0;
    for (int i = 0; i < N; ++i, y *= x) {
        z[i] = y;
    }
    return z;
}

template <int N>
typename Vandermonde<N>::RowVector
Vandermonde<N>::differentiate1(double x) {
    RowVector z = RowVector::Zero();
    double y = 1.0;
    for (int i = 1; i < N; ++i, y *= x) {
        z[i] = i*y;
    }
    return z;
}

template <int N>
typename Vandermonde<N>::RowVector
Vandermonde<N>::differentiate2(double x) {
    RowVector z = RowVector::Zero();
    double y = 1.0;
    for (int i = 2; i < N; ++i, y *= x) {
        z[i] = i*(i-1)*y;
    }
    return z;
}

template <int N>
typename Vandermonde<N>::RowVector
Vandermonde<N>::moment(double x0, double x1, int m) {
    RowVector z = RowVector::Zero();
    double y0 = x0;
    double y1 = x1;
    for (int j = 0; j < m; ++j, y0 *= x0, y1 *= x1);
    for (int i = 0; i < N; ++i, y0 *= x0, y1 *= x1) {
        z[i] = (y1 / (i+m+1)) - (y0 / (i+m+1));
    }
    return z;
}

Eigen::Vector4d solveRampPoly(double v0, double v1, double x0, double x1, double s0, double s1) {
    Eigen::Vector4d b;
    Eigen::Matrix4d m;
    // p(x0) = v0
    m.row(0) = Vandermonde<4>::eval(x0);
    b[0] = v0;
    // p(x1) = v1
    m.row(1) = Vandermonde<4>::eval(x1);
    b[1] = v1;
    // p'(x0) = s0
    m.row(2) = Vandermonde<4>::differentiate1(x0);
    b[2] = s0;
    // p'(x1) = s1
    m.row(3) = Vandermonde<4>::differentiate1(x1);
    b[3] = s1;
    return m.fullPivLu().solve(b);
}

template class Vandermonde<4>;

}}}} // namespace lsst::meas::modelfit::detail
