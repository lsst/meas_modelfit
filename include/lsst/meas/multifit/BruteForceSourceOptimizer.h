// -*- LSST-C++ -*-
/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010, 2011 LSST Corporation.
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

#ifndef LSST_MEAS_MULTIFIT_BruteForceSourceOptimizer
#define LSST_MEAS_MULTIFIT_BruteForceSourceOptimizer

#include "lsst/meas/multifit/GaussianDistribution.h"
#include "lsst/meas/multifit/Evaluator.h"

namespace lsst {
namespace meas {
namespace multifit {

class BruteForceSourceOptimizer {
public:

    BruteForceSourceOptimizer() {}

    void solve(Evaluator::Ptr const & evaluator, int n);

    int getBestIndex() const { return _bestIndex; }
    double getBestObjectiveValue() const { return _objectiveValues[_bestIndex]; }
    lsst::ndarray::Array<double const,1,1> getBestParameters() const { return _parameters[_bestIndex]; }
    lsst::ndarray::Array<double const,1,1> getBestCoefficients() const { return _bestCoefficients; }
    lsst::ndarray::Array<double const,2,2> getCoefficientCovariance() const { return _coefficientCovariance; }

    lsst::ndarray::Array<double const,2,2> getParameters() const { return _parameters; }
    lsst::ndarray::Array<double const,1,1> getObjectiveValues() const { return _objectiveValues; }

private:
    int _bestIndex;
    lsst::ndarray::Array<double,2,2> _parameters;
    lsst::ndarray::Array<double,1,1> _objectiveValues;
    lsst::ndarray::Array<double,2,2> _coefficientCovariance;
    lsst::ndarray::Array<double,1,1> _bestCoefficients;
};

}}} //end namespace lsst::meas::multifit

#endif //end #ifndef LSST_MEAS_MULTIFIT_BruteForceSourceOptimizer
