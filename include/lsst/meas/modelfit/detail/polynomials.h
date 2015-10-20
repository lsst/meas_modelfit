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

#ifndef LSST_MEAS_MODELFIT_DETAIL_polynomials_h_INCLUDED
#define LSST_MEAS_MODELFIT_DETAIL_polynomials_h_INCLUDED

#include "Eigen/Core"

namespace lsst { namespace meas { namespace modelfit { namespace detail {

/**
 * Class that computes rows of the Vandermonde matrix and related matrices;
 * the dot product of these row vectors with the polynomial coefficient
 * vectors evaluates the polynomial (or computes a derivative).
*/
template <int N>
struct Vandermonde {

    typedef Eigen::Matrix<double,1,N> RowVector;

    /// Return a row vector that product with a polynomial coefficient vector[
    /// evaluates the polynomial at x.
    static RowVector eval(double x);

    /// Return a row vector whose product with a polynomial coefficient vector
    /// evaluates the first derivative at x.
    static RowVector differentiate1(double x);

    /// Return a row vector whose product with a polynomial coefficient vector
    /// evaluates the second derivative at x.
    static RowVector differentiate2(double x);

    /// Return a row vector whose product with a polynomial coefficient vector
    /// computes the integral of p(x) x^m dx from x0 to x1
    static RowVector moment(double x0, double x1, int m=0);

};

/// Solve for the coefficients of a cubic polynomial p(x) that goes from
/// p(x0)=v0 to p(x1)=v1, with p'(x0)=s0 and p'(x1)=s1.
Eigen::Vector4d solveRampPoly(double v0, double v1, double x0, double x1, double s0, double s1);

}}}} // namespace lsst::meas::modelfit::detail

#endif // !LSST_MEAS_MODELFIT_DETAIL_polynomials_h_INCLUDED
