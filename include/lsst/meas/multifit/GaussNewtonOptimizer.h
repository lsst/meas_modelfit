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

#ifndef LSST_MEAS_MULTIFIT_GaussNewtonOptimizer
#define LSST_MEAS_MULTIFIT_GaussNewtonOptimizer

#include "lsst/meas/multifit/BaseEvaluator.h"
#include "lsst/meas/multifit/SimpleInterpreter.h"
#include "lsst/meas/multifit/GaussianDistribution.h"

namespace lsst {
namespace meas {
namespace multifit {

class GaussNewtonOptimizer {
public:
    GaussNewtonOptimizer(){}

    GaussianDistribution::Ptr solve(
        BaseEvaluator::ConstPtr const & evaluator,
        double const fTol=1.e-8, double const gTol=1.e-8, 
        double const minStep=1.e-8, 
        int const maxIter=200, 
        double const tau=1.e-3, 
        bool retryWithSvd=false
    );

    
    ndarray::Array<const double, 1, 1> getParameters() const;
    ndarray::Array<const double, 1, 1> getCoefficients() const
    

private:

    Evaluation _evaluation;
    bool _solverSuccess;
};

}}} //end namespace lsst::meas::multifit

};


#endif //end #ifndef LSST_MEAS_MULTIFIT_UnifiedNLSolver
