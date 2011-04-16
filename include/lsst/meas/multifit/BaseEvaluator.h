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

#ifndef LSST_MEAS_MULTIFIT_BaseEvaluator
#define LSST_MEAS_MULTIFIT_BaseEvaluator

#include "lsst/ndarray.h"
#include "lsst/meas/multifit/constants.h"
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

namespace lsst { namespace meas { namespace multifit {

/**
 *  @brief An abstract base class that combines a data vector and model(s),
 *         organizing the model as a parameterized matrix multiplied by a vector
 *         of linear coefficients.
 *
 *  BaseEvaluator subclasses should be immutable.
 */
class BaseEvaluator {
public:

    typedef boost::shared_ptr<BaseEvaluator> Ptr;

    /// @brief Size of data vector (number of rows of matrix).
    int getDataSize() const { return _dataVector.getSize<0>(); }

    /// @brief Size of coefficient vector (number of colums of matrix).
    int getCoefficientSize() const { return _coefficientSize; }

    /// @brief Number of parameters.
    int getParameterSize() const { return _parameterSize; }

    /**
     *  @brief Return the sum of the log of the variance of all data points.
     *
     *  More precisely, if the variance for pixel @f$i@f$ is @f$\sigma_i^2@f$, this function
     *  should return @f$ \sum_i \ln \sigma_i^2 @f$.
     *
     *  This can be used to transform the @f$\chi^2@f$ into the log likelihood, which can
     *  be useful for certain Bayesian applications.
     */
    double getLogVarianceSum() const { return _logVarianceSum; }

    /**
     *  @brief Data vector.
     *
     *  If the data vector is weighted (divided by sigma) the evaluted model matrix should be as well.
     */
    lsst::ndarray::Array<Pixel const,1,1> getDataVector() const { return _dataVector; }

    /**
     *  @brief Evaluate the matrix with the given parameters.
     *
     *  @param[out] matrix  An array to fill with shape (getDataSize(), getCoefficientSize()).
     *  @param[in]  param   An array of parameters with size getParameterSize().
     *
     *  If the data vector is weighted, the output matrix should be as well (each row should be divided
     *  by the corresponding pixel sigma value).
     */
    void evaluateModelMatrix(
        lsst::ndarray::Array<Pixel,2,2> const & matrix,
        lsst::ndarray::Array<Pixel const,1,1> const & param
    ) const;

    /**
     *  @brief Evaluate the derivative of the matrix with respect to the parameters.
     *
     *  The first dimension of the array is the parameter dimension.
     */
    void evaluateModelDerivative(
        lsst::ndarray::Array<Pixel,3,3> const & derivative,
        lsst::ndarray::Array<Pixel const,1,1> const & param
    ) const;

    void writeInitialParameters(lsst::ndarray::Array<Pixel,1,1> const & param) const;

    virtual ~BaseEvaluator() {}

protected:

    BaseEvaluator(int dataSize, int coefficientSize, int parameterSize, double logVarianceSum) :
        _coefficientSize(coefficientSize),
        _parameterSize(parameterSize),
        _dataVector(ndarray::allocate(ndarray::makeVector(dataSize))),
        _logVarianceSum(logVarianceSum)
    {}

    BaseEvaluator(
        ndarray::Array<Pixel,1,1> const & data, int coefficientSize, int parameterSize, 
        double logVarianceSum
    ) :
        _coefficientSize(coefficientSize),
        _parameterSize(parameterSize),
        _dataVector(data),
        _logVarianceSum(logVarianceSum)
    {}

    BaseEvaluator(BaseEvaluator const & other) :
        _coefficientSize(other._coefficientSize), 
        _parameterSize(other._parameterSize), 
        _dataVector(other._dataVector),
        _logVarianceSum(other._logVarianceSum)
    {}

    virtual void _evaluateModelMatrix(
        ndarray::Array<Pixel,2,2> const & matrix,
        ndarray::Array<Pixel const,1,1> const & param
    ) const = 0;

    virtual void _evaluateModelDerivative(
        ndarray::Array<Pixel,3,3> const & derivative,
        ndarray::Array<Pixel const,1,1> const & param
    ) const;

    virtual void _writeInitialParameters(ndarray::Array<Pixel,1,1> const & param) const = 0;

    int const _coefficientSize;
    int const _parameterSize;
    ndarray::Array<Pixel,1,1> _dataVector;
    double _logVarianceSum;

private:
    void operator=(BaseEvaluator const &) {}
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_BaseEvaluator
