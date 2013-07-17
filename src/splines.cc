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

#include "Eigen/QR"
#include "boost/noncopyable.hpp"
#include "boost/format.hpp"
#include "gsl/gsl_bspline.h"
#include "gsl/gsl_vector.h"

#include "ndarray/eigen.h"
#include "lsst/pex/exceptions.h"
#include "lsst/meas/multifit/splines.h"

namespace lsst { namespace meas { namespace multifit {

namespace {

void checkStatus(int status, std::string const & message) {
    switch (status) {
    case 0:
        return;
    case GSL_EDOM:
        throw LSST_EXCEPT(
            pex::exceptions::DomainErrorException,
            message + "GSL domain error"
        );
    case GSL_ERANGE:
        throw LSST_EXCEPT(
            pex::exceptions::OverflowErrorException,
            message + "GSL overflow/underflow error"
        );
    case GSL_ENOMEM:
        throw LSST_EXCEPT(
            pex::exceptions::MemoryException,
            message + "GSL memory allocation error"
        );
    case GSL_EINVAL:
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            message + "GSL invalid argument error"
        );
    default:
        throw LSST_EXCEPT(
            pex::exceptions::RuntimeErrorException,
            message + "unknown GSL error"
        );
    }
}

} // anonymous

class SplineBasis::Impl : private boost::noncopyable {
public:

    explicit Impl(int nKnots, int degree) {
        _w = gsl_bspline_alloc(degree + 1, nKnots);
        _dw = gsl_bspline_deriv_alloc(degree + 1);
    }

    ~Impl() {
        gsl_bspline_free(_w);
        gsl_bspline_deriv_free(_dw);
    }

    double _min;
    double _max;
    gsl_bspline_workspace * _w;
    gsl_bspline_deriv_workspace * _dw;
};

SplineBasis::SplineBasis(ndarray::Array<double const,1,1> const & knots, int degree) :
    _impl(new Impl(knots.getSize<0>(), degree))
{
    if (knots.getSize<0>() < 2) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            "Knots vector must have at least 2 elements"
        );
    }
    for (int i = 1, n = knots.getSize<0>(); i < n; ++i) {
        if (knots[i] < knots[i - 1]) {
            throw LSST_EXCEPT(
                pex::exceptions::LogicErrorException,
                "Knots vector must be nondecreasing"
            );
        }
    }
    _impl->_min = knots[0];
    _impl->_max = knots[knots.getSize<0>() - 1];
    gsl_vector_const_view gknots = gsl_vector_const_view_array(knots.getData(), knots.getSize<0>());
    gsl_bspline_knots(&gknots.vector, _impl->_w);
}

int SplineBasis::getBasisSize() const {
    return gsl_bspline_ncoeffs(_impl->_w);
}

int SplineBasis::getOrder() const {
    return gsl_bspline_order(_impl->_w);
}

ndarray::Array<double const,1,0> SplineBasis::getKnots() const {
    ndarray::Array<double const,1,0> knots = ndarray::external(
        _impl->_w->knots->data,
        ndarray::Vector<int,1>(_impl->_w->knots->size),
        ndarray::Vector<int,1>(_impl->_w->knots->stride),
        _impl
    );
    return knots;
}

void SplineBasis::evaluate(double x, ndarray::Array<double,1,1> const & b) const {
    if (b.getSize<0>() != getBasisSize()) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthErrorException,
            (boost::format("Basis size (%d) does not match columns of b (%d)")
             % getBasisSize() % b.getSize<1>()).str()
        );
    }
    if (x < _impl->_min || x > _impl->_max) {
        b.deep() = 0.0;
    } else {
        gsl_vector_view gb = gsl_vector_view_array(b.getData(), b.getSize<0>());
        checkStatus(gsl_bspline_eval(x, &gb.vector, _impl->_w), "In basis spline evaluation: ");
    }
}

ndarray::Array<double,1,1> SplineBasis::evaluate(double x) const {
    ndarray::Array<double,1,1> b = ndarray::allocate(getBasisSize());
    evaluate(x, b);
    return b;
}

void SplineBasis::evaluate(
    ndarray::Array<double const,1,1> const & x,
    ndarray::Array<double,2,1> const & b
) const {
    if (x.getSize<0>() != b.getSize<0>()) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthErrorException,
            (boost::format("Size of x (%d) does not match rows of b (%d)")
             % x.getSize<0>() % b.getSize<0>()).str()
        );
    }
    if (b.getSize<1>() != getBasisSize()) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthErrorException,
            (boost::format("Basis size (%d) does not match columns of b (%d)")
             % getBasisSize() % b.getSize<1>()).str()
        );
    }
    gsl_matrix_view gb = gsl_matrix_view_array_with_tda(
        b.getData(), b.getSize<0>(), b.getSize<1>(), b.getStride<0>()
    );
    int const nx = x.getSize<0>();
    for (int i = 0; i < nx; ++i) {
        if (x[i] < _impl->_min || x[i] > _impl->_max) {
            b[i] = 0.0;
        } else {
            gsl_vector_view gbi = gsl_matrix_row(&gb.matrix, i);
            checkStatus(gsl_bspline_eval(x[i], &gbi.vector, _impl->_w), "In basis spline evaluation: ");
        }
    }
}

ndarray::Array<double,2,2> SplineBasis::evaluate(
    ndarray::Array<double const,1,1> const & x
) const {
    ndarray::Array<double,2,2> b = ndarray::allocate(x.getSize<0>(), getBasisSize());
    evaluate(x, b);
    return b;
}


void SplineBasis::evaluateDerivatives(double x, ndarray::Array<double,2,1> const & db) const {
    if (db.getSize<0>() != getBasisSize()) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthErrorException,
            (boost::format("Basis size (%d) does not match rows of db (%d)")
             % getBasisSize() % db.getSize<0>()).str()
        );
    }
    gsl_matrix_view gdb = gsl_matrix_view_array_with_tda(
        db.getData(), db.getSize<0>(), db.getSize<1>(), db.getStride<0>()
    );
    checkStatus(gsl_bspline_deriv_eval(x, db.getSize<1>()-1, &gdb.matrix, _impl->_w, _impl->_dw),
                "In basis spline derivative evaluation: ");
}

ndarray::Array<double,2,2> SplineBasis::evaluateDerivatives(double x, int nDeriv) const {
    ndarray::Array<double,2,2> db = ndarray::allocate(getBasisSize(), nDeriv + 1);
    evaluateDerivatives(x, db);
    return db;
}

void SplineBasis::evaluateIntegral(ndarray::Array<double,1,1> const & ib) const {
    if (ib.getSize<0>() != getBasisSize()) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthErrorException,
            (boost::format("Basis size (%d) does not match columns of ib (%d)")
             % getBasisSize() % ib.getSize<1>()).str()
        );
    }
    int const order = getOrder(), basisSize = getBasisSize();
    for (int i = 0; i < basisSize; ++i) {
        ib[i] = (gsl_vector_get(_impl->_w->knots, i + order) - gsl_vector_get(_impl->_w->knots, i)) / order;
    }
}

ndarray::Array<double,1,1> SplineBasis::evaluateIntegral() const {
    ndarray::Array<double,1,1> ib = ndarray::allocate(getBasisSize());
    evaluateIntegral(ib);
    return ib;
}

SplineBasis::~SplineBasis() {}



ConstrainedSplineBasis::ConstrainedSplineBasis(SplineBasis const & spline) :
    _spline(spline), _constraintMatrix(), _constraintVector(),
    _y1(Eigen::VectorXd::Zero(spline.getBasisSize())),
    _q2(Eigen::MatrixXd::Identity(spline.getBasisSize(), spline.getBasisSize()))
{}


void ConstrainedSplineBasis::evaluate(double x, ndarray::Array<double,1,1> const & b, double & c) const {
    ndarray::Array<double,1,1> sb = _spline.evaluate(x);
    b.asEigen() = sb.asEigen() * _q2;
    c = sb.asEigen().dot(_y1);
}

void ConstrainedSplineBasis::evaluate(
    ndarray::Array<double const,1,1> const & x,
    ndarray::Array<double,2,1> const & b,
    ndarray::Array<double,1,1> const & c
) const {
    ndarray::Array<double,2,2> sb = _spline.evaluate(x);
    b.asEigen() = sb.asEigen() * _q2;
    c.asEigen() = sb.asEigen() * _y1;
}

void ConstrainedSplineBasis::unconstrainCoefficients(
    ndarray::Array<double const,1,1> const & constrained,
    ndarray::Array<double,1,1> const & unconstrained
) const {
    if (constrained.getSize<0>() != getBasisSize()) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthErrorException,
            (boost::format("Constrained coefficient vector size (%d) does not match basis size (%d)")
             % constrained.getSize<0>() % getBasisSize()).str()
        );
    }
    if (unconstrained.getSize<0>() != _spline.getBasisSize()) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthErrorException,
            (boost::format("Unconstrained coefficient vector size (%d) does not match spline basis size (%d)")
             % unconstrained.getSize<0>() % _spline.getBasisSize()).str()
        );
    }
    unconstrained.asEigen() = _y1 + _q2 * constrained.asEigen();
}

ndarray::Array<double,1,1> ConstrainedSplineBasis::unconstrainCoefficients(
    ndarray::Array<double const,1,1> const & constrained
) const {
    ndarray::Array<double,1,1> unconstrained = ndarray::allocate(_spline.getBasisSize());
    unconstrainCoefficients(constrained, unconstrained);
    return unconstrained;
}

void ConstrainedSplineBasis::addConstraint(double x, double v, int n) {
    ndarray::Array<double,2,-2> derivatives = _spline.evaluateDerivatives(x, n).transpose();
    addManualConstraint(derivatives[n], v);
}

void ConstrainedSplineBasis::addIntegralConstraint(double i) {
    ndarray::Array<double,1,1> k = _spline.evaluateIntegral();
    addManualConstraint(k, i);
}

void ConstrainedSplineBasis::addManualConstraint(ndarray::Array<double const,1,0> const & k, double c) {
    int nc = _constraintMatrix.getSize<0>();
    int ns = _spline.getBasisSize();
    if (k.getSize<0>() != ns) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthErrorException,
            (boost::format("Size of constraint vector (%d) does not match unconstrained basis size (%d)")
             % k.getSize<0>() % ns).str()
        );
    }
    ndarray::Array<double,2,2> newConstraintMatrix = ndarray::allocate(nc + 1, ns);
    ndarray::Array<double,1,1> newConstraintVector = ndarray::allocate(nc + 1);
    if (nc > 0) {
        newConstraintMatrix[ndarray::view(0, nc)] = _constraintMatrix;
        newConstraintVector[ndarray::view(0, nc)] = _constraintVector;
    }
    double scaling = std::max(k.asEigen().lpNorm<Eigen::Infinity>(), std::abs(c));
    newConstraintMatrix[nc] = k / scaling;
    newConstraintVector[nc] = c / scaling;
    ++nc;
    newConstraintMatrix.swap(_constraintMatrix);
    newConstraintVector.swap(_constraintVector);
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(_constraintMatrix.asEigen().adjoint());
    Eigen::MatrixXd q = qr.householderQ();
    Eigen::MatrixXd r = qr.matrixQR().triangularView<Eigen::Upper>();
    assert(r.cols() == nc);
    assert(r.rows() == ns);
    Eigen::VectorXd z1 = _constraintVector.asEigen();
    r.topRows(nc).triangularView<Eigen::Upper>().adjoint().solveInPlace(z1);
    _y1 = q.leftCols(nc) * z1;
    _q2 = q.rightCols(ns - nc);
}

double SplineFunction::operator()(double x) const {
    _basis.evaluate(x, _ws);
    return _ws.asEigen().dot(_coefficients.asEigen());
}

void SplineFunction::operator()(
    ndarray::Array<double const,1,0> const & x, ndarray::Array<double,1,0> const & v
) const {
    if (x.getSize<0>() != v.getSize<0>()) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthErrorException,
            (boost::format("Sizes of x (%d) and v (%d) do not agree")
             % x.getSize<0>() % v.getSize<0>()).str()
        );
    }
    ndarray::Array<double const,1,0>::Iterator ix = x.begin(), xend = x.end();
    ndarray::Array<double,1,0>::Iterator iv = v.begin();
    for (; ix != xend; ++ix, ++iv) {
        *iv = (*this)(*ix);
    }
}

ndarray::Array<double,1,1> SplineFunction::operator()(ndarray::Array<double const,1,1> const & x) const {
    ndarray::Array<double,1,1> v = ndarray::allocate(x.getShape());
    (*this)(x, v);
    return v;
}



}}} // namespace lsst::meas::multifit
