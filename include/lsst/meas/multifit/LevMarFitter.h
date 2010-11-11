// -*- lsst-c++ -*-

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
 
/**
 * @file
 *
 * Declaration levmar-based fitter.
 */
#ifndef LSST_MEAS_MULTIFIT_LEVMARFITTER_H
#define LSST_MEAS_MULTIFIT_LEVMARFITTER_H

#include <float.h>
#include "lsst/pex/policy/DefaultPolicyFile.h"
#include "lsst/pex/policy/Policy.h"
#include "lsst/pex/exceptions/Runtime.h"
#include "lsst/meas/multifit/ModelEvaluator.h"

namespace lsst {
namespace meas {
namespace multifit {

struct LevMarFitterResult {
public:

    enum TerminationEnum {
        NONE = 0,
        GRADIENT_SMALL = 1,
        PARAMETERS_UNCHANGED = 2,
        MAX_ITERATION_REACHED = 3,
        SINGULAR_MATRIX = 4,
        NO_REDUCTION_POSSIBLE = 5,
        RESIDUALS_SMALL = 6,
        INVALID_FUNCTION_VALUES = 7
    };

#ifndef SWIG
    LevMarFitterResult(
        double * levmarInfo,
        Model::ConstPtr const & model,
        Eigen::VectorXd const & parameters,
        Eigen::MatrixXd const & covariance
    );
#endif
    
    LevMarFitterResult() : 
        termination(NONE), chisqInitial(DBL_MAX), chisqFinal(DBL_MAX),
        nIterations(0), nFunctionEvaluations(0), nJacobianEvaluations(0) {}

    bool hasConverged() const {
        return termination == GRADIENT_SMALL || termination == PARAMETERS_UNCHANGED
            || termination == RESIDUALS_SMALL;
    }

    bool needsRestart() const {
        return termination == SINGULAR_MATRIX || termination == NO_REDUCTION_POSSIBLE;
    }

    TerminationEnum termination;
    double chisqInitial;
    double chisqFinal;
    double maxGradient;
    double lastStepNorm;
    int nIterations;
    int nFunctionEvaluations;
    int nJacobianEvaluations;
    int nMatrixFactorizations;
    Model::ConstPtr model;
    Eigen::VectorXd parameters;
    Eigen::MatrixXd covariance;
};

class LevMarFitter {
public:

    typedef LevMarFitterResult Result;

    LevMarFitter(
        lsst::pex::policy::Policy::Ptr const & policy = lsst::pex::policy::Policy::Ptr()
    );
    ~LevMarFitter() {} 
    Result apply(multifit::ModelEvaluator::Ptr const & evaluator) const;

    Eigen::VectorXd checkDerivatives(multifit::ModelEvaluator::Ptr const & evaluator) const;

    /**
     * Retrieve the default policy for configuring a LevMarFitter
     *
     * The defaults are provided in the dictionary file
     */
    static PTR(lsst::pex::policy::PolicySource) getDefaultPolicySource() {
        static const PTR(lsst::pex::policy::PolicySource) source(
            new lsst::pex::policy::DefaultPolicyFile(
                "meas_multifit", 
                "LevMarFitterDict.paf", 
                "policy"
            )
        );
        return source;
    }

private:
    lsst::pex::policy::Policy::Ptr _policy;
};

}}}

#endif
