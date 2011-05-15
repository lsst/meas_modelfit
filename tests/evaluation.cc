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
#define BOOST_TEST_MODULE evaluation

#include "boost/test/unit_test.hpp"
#include "lsst/meas/multifit/constants.h"
#include "lsst/meas/multifit/Evaluation.h"
#include "lsst/ndarray/eigen.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/detection.h"
#include <Eigen/Array>


namespace lsst { namespace meas { namespace multifit {

class TestEvaluator : public BaseEvaluator {
public:

    virtual double clipToBounds(lsst::ndarray::Array<double,1,1> const & parameters) const { return 0.0; }

    TestEvaluator(int dataSize, int coefficientSize) :
        BaseEvaluator(dataSize, coefficientSize, 0),
        _matrix(ndarray::allocate(dataSize, coefficientSize))
    {
        ndarray::viewAsEigen(_matrix) = Eigen::MatrixXd::Random(dataSize, coefficientSize);
        ndarray::viewAsEigen(_dataVector) = Eigen::VectorXd::Random(dataSize);
    }

protected:
    
    virtual void _writeInitialParameters(ndarray::Array<Pixel,1,1> const & param) const {}

    virtual void _evaluateModelMatrix(
        ndarray::Array<Pixel,2,2> const & matrix,
        ndarray::Array<Pixel const,1,1> const & param
    ) const {
        matrix.deep() = _matrix;
    }

    ndarray::Array<Pixel,2,2> _matrix;
};

}}}

using namespace lsst::meas::multifit;

BOOST_AUTO_TEST_CASE(evaluation1) {
    for (int n = 0; n < 1; ++n) {
        TestEvaluator::Ptr evaluator(new TestEvaluator(50, 10));
        Evaluation robustEv(evaluator);
        robustEv.getCoefficients();
    }
}
