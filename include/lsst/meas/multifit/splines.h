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

#ifndef LSST_MEAS_MULTIFIT_splines_h_INCLUDED
#define LSST_MEAS_MULTIFIT_splines_h_INCLUDED

#include "lsst/base.h"
#include "ndarray.h"

namespace lsst { namespace meas { namespace multifit {

class SplineBasis {
public:

    explicit SplineBasis(ndarray::Array<double const,1,1> const & knots, int degree=3);

    /// Return the number of basis functions
    int getBasisSize() const;

    /// Return the order of the spline (one greater than the degree; cubic spline has order==4).
    int getOrder() const;

    /// Return the degree of the spline (one less than the order; cubic spline has degree==3).
    int getDegree() const { return getOrder() - 1; }

    /// Return the vector of knots (includes repeated knots at the beginning and end).
    ndarray::Array<double const,1,0> getKnots() const;

    /// Evaluate the spline basis at the given point, into the given vector.
    void evaluate(double x, ndarray::Array<double,1,1> const & b) const;

    /// Evaluate the spline basis at the given point.
    ndarray::Array<double,1,1> evaluate(double x) const;

    /// @brief Evaluate the spline basis at the given points into a (x.size(), getBasisSize()) matrix.
    void evaluate(
        ndarray::Array<double const,1,1> const & x,
        ndarray::Array<double,2,1> const & b
    ) const;

    /// @brief Evaluate the spline basis at the given points into a (x.size(), getBasisSize()) matrix.
    ndarray::Array<double,2,2> evaluate(
        ndarray::Array<double const,1,1> const & x
    ) const;

    /**
     *  @brief Evaluate the derivatives of the spline basis functions at the given point, into a
     *         (getBasisSize(), nDeriv+1) matrix.
     *
     *  The first column of the matrix is always the function itself and intermediate derivatives
     *  are always computed in subsequent columns, so to evaluate e.g. the second derivative,
     *  a matrix with 3 columns should be provided, and then the last column inspected.
     */
    void evaluateDerivatives(double x, ndarray::Array<double,2,1> const & db) const;

    /**
     *  @brief Evaluate the derivatives of the spline basis functions at the given point, into a
     *         (getBasisSize(), nDeriv+1) matrix.
     *
     *  The first column of the matrix is always the function itself and intermediate derivatives
     *  are always computed, so if nDeriv=2, a matrix with 3 columns will be be returned with the
     *  second derivative in the last column.
     */
    ndarray::Array<double,2,2> evaluateDerivatives(double x, int nDeriv) const;

    /// @brief Evaluate the definite integral of the basis functions over all real numbers.
    void evaluateIntegral(ndarray::Array<double,1,1> const & ib) const;

    /// @brief Evaluate the definite integral of the basis functions over all real numbers.
    ndarray::Array<double,1,1> evaluateIntegral() const;

    ~SplineBasis(); // needs to be in .cc so compiler can see ~Impl()

private:
    class Impl;
    PTR(Impl) _impl;
};

class ConstrainedSplineBasis {
public:

    explicit ConstrainedSplineBasis(SplineBasis const & spline);

    SplineBasis const & getSplineBasis() const { return _spline; }

    int getBasisSize() const { return _q2.cols(); }

    void evaluate(double x, ndarray::Array<double,1,1> const & b, double & c) const;

    void evaluate(
        ndarray::Array<double const,1,1> const & x,
        ndarray::Array<double,2,1> const & b,
        ndarray::Array<double,1,1> const & c
    ) const;

    void unconstrainCoefficients(
        ndarray::Array<double const,1,1> const & constrained,
        ndarray::Array<double,1,1> const & unconstrained
    ) const;

    ndarray::Array<double,1,1> unconstrainCoefficients(
        ndarray::Array<double const,1,1> const & constrained
    ) const;

    /**
     *  @brief Add a constraint that fixes the nth derivative at x to v.
     *
     *  When n is zero, fixes the value of the spline basis at x to v.
     */
    void addConstraint(double x, double v=0.0, int n=0);

    /**
     *  @brief Add a constraint that fixes the full integral to the given value.
     */
    void addIntegralConstraint(double i=1.0);

    /**
     *  @brief Add a constraint that fixes the dot product of k with the coefficients to v.
     *
     *  The size of k must be the same as the size of the original (unconstrained) spline basis.
     */
    void addManualConstraint(ndarray::Array<double const,1,0> const & k, double v);

private:
    SplineBasis _spline;
    ndarray::Array<double,2,2> _constraintMatrix;
    ndarray::Array<double,1,1> _constraintVector;
    Eigen::VectorXd _y1;
    Eigen::MatrixXd _q2;
};

class SplineFunction {
public:

    SplineFunction(SplineBasis const & basis, ndarray::Array<double const,1,1> const & coefficients) :
        _basis(basis), _coefficients(coefficients), _ws(ndarray::allocate(coefficients.getShape()))
    {}

    double operator()(double x) const;

    void operator()(ndarray::Array<double const,1,0> const & x, ndarray::Array<double,1,0> const & v) const;

    ndarray::Array<double,1,1> operator()(ndarray::Array<double const,1,1> const & x) const;

private:
    SplineBasis _basis;
    ndarray::Array<double const,1,1> _coefficients;
    ndarray::Array<double,1,1> _ws;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_splines_h_INCLUDED
