// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
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
 
#include <iostream>
#include <cmath>
#include <vector>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE qp

#include "boost/test/unit_test.hpp"
#include "boost/test/floating_point_comparison.hpp"
#include "lsst/meas/multifit/qp.h"
#include "lsst/ndarray/eigen.h"

#include <Eigen/Array>

namespace nd = lsst::ndarray;

namespace {

nd::Array<double,2,1> makeRandomMatrix(int rows, int cols) {
    nd::Array<double,2,2> r = nd::allocate(rows, cols);
    if (rows * cols > 0)
        nd::viewAsEigen(r).setRandom();
    return r;
}

nd::Array<double,1,1> makeRandomVector(int rows) {
    nd::Array<double,1,1> r = nd::allocate(rows);
    if (rows > 0)
        nd::viewAsEigen(r).setRandom();
    return r;
}

double checkQP(
    nd::Array<double,2> const & g,
    nd::Array<double,1> const & c,
    nd::Array<double,2> const & ae,
    nd::Array<double,1> const & be,
    nd::Array<double,2> const & ai,
    nd::Array<double,1> const & bi,
    nd::Array<double,1> const & x
) {
    Eigen::IOFormat fmt(16);
    double r = 0.5 * nd::viewAsEigen(x).dot(
        nd::viewAsEigen(g) * nd::viewAsEigen(x)
    ) + nd::viewAsEigen(c).dot(nd::viewAsEigen(x));
    if (be.size() > 0) {
        if (!nd::viewAsEigen(be).isApprox(nd::viewAsEigen(ae) * nd::viewAsEigen(x))) {
            r = std::numeric_limits<double>::infinity();
        }
    }
    if (bi.size() > 0) {
        if (!((nd::viewAsEigen(ai) * nd::viewAsEigen(x) - nd::viewAsEigen(bi)).cwise() >= 
              -std::numeric_limits<double>::epsilon()).all()) {
            r = std::numeric_limits<double>::infinity();
        }
    }
    return r;
}

void testQP(
    int const nx,
    int const nd,
    int const ni,
    int const ne
) {
    int const MAX_ITER = 100;
    double const SQRT_EPS = std::sqrt(std::numeric_limits<double>::epsilon());

    bool success = false;
    int iterations = 0;

    while (!success && ++iterations < MAX_ITER) {
        nd::Array<double,2> g = nd::allocate(nx, nx);
        Eigen::MatrixXd j = Eigen::MatrixXd::Random(nd, nx);
        nd::viewAsEigen(g).part<Eigen::SelfAdjoint>() = j.transpose() * j;

        nd::Array<double,1> c = makeRandomVector(nx);

        nd::Array<double,2> ai = makeRandomMatrix(ni, nx);
        nd::Array<double,1> bi = makeRandomVector(ni);

        nd::Array<double,2> ae = makeRandomMatrix(ne, nx);
        nd::Array<double,1> be = makeRandomVector(ne);

        nd::Array<double,1> x = nd::allocate(nx);

        double r = lsst::meas::multifit::QPSolver(g, c).equality(ae, be).inequality(ai, bi).solve(x);
        if (r < std::numeric_limits<double>::infinity()) {
            double s = checkQP(g, c, ae, be, ai, bi, x);
            BOOST_CHECK_CLOSE(s, r, 1E-8);
            for (int i = 0; i < nx; ++i) {
                x[i] += SQRT_EPS;
                BOOST_CHECK( checkQP(g, c, ae, be, ai, bi, x) > r );
                x[i] -= 2.0 * SQRT_EPS;
                BOOST_CHECK( checkQP(g, c, ae, be, ai, bi, x) > r );
                x[i] += SQRT_EPS;
            }

            success = true;
        }
    }
}

} // anonymous

BOOST_AUTO_TEST_CASE(qp) {
    testQP(12, 20, 4, 3);
    testQP(12, 20, 0, 3);
    testQP(12, 20, 4, 0);
}
