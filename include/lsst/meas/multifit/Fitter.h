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

#ifndef LSST_MEAS_MULTIFIT_Fitter_h_INCLUDED
#define LSST_MEAS_MULTIFIT_Fitter_h_INCLUDED

#include "lsst/meas/multifit/constants.h"

namespace lsst { namespace meas { namespace multifit {

class JointObjective {
public:

    JointObjective(int dataSize, int parameterSize) :
        _dataSize(dataSize), _parameterSize(parameterSize)
    {}

    int getDataSize() const { return _dataSize; }
    int getParameterSize() const { return _parameterSize; }

    virtual void computeResiduals(
        ndarray::Array<double const,1,1> const & parameters,
        ndarray::Array<double,1,1> const & residuals
    ) const = 0;

    virtual bool hasPrior() const { return false; }

    virtual double computePrior(ndarray::Array<double const,1,1> const & parameters) const;

    virtual void differentiatePrior(
        ndarray::Array<double const,1,1> const & parameters,
        ndarray::Array<double,1,1> const & gradient,
        ndarray::Array<double,2,1> const & hessian
    ) const;

    virtual ~JointObjective() {}

protected:
    int _dataSize;
    int _parameterSize;
};

class MarginalObjective {
public:

    MarginalObjective(int dataSize, int parameterSize, int nestedSize) :
        _dataSize(dataSize), _parameterSize(parameterSize), _nestedSize(nestedSize)
    {}

    int getDataSize() const { return _dataSize; }
    int getParameterSize() const { return _parameterSize; }
    int getNestedSize() const { return _nestedSize; }

    virtual double getDataSquaredNorm() const = 0;

    virtual LogGaussian evaluate(
        ndarray::Array<double const,1,1> const & parameters
    ) const = 0;

    virtual ~MarginalObjective() {}

};

class Fitter {
public:

    

    virtual ~Fitter() {}
};

class FitResult {
public:

    virtual ~FitResult() {}
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_constants_h_INCLUDED
