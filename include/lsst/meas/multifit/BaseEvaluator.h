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

#ifndef LSST_MEAS_MULTIFIT_BaseEvaluator
#define LSST_MEAS_MULTIFIT_BaseEvaluator

#include "lsst/ndarray.hpp"
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
    int getDataSize() const { return _data_vector.getSize<0>(); }

    /// @brief Size of coefficient vector (number of colums of matrix).
    int getCoefficientSize() const { return _coefficient_size; }

    /// @brief Number of parameters.
    int getParameterSize() const { return _parameter_size; }

    /// @brief Data vector.
    ndarray::Array<double const,1,1> getDataVector() const { return _data_vector; }

    /// @brief Evaluate the matrix with the given parameters.
    void evaluateModelMatrix(
        ndarray::Array<double,2,2> const & matrix,
        ndarray::Array<double const,1,1> const & param
    ) const;

    /**
     *  @brief Evaluate the derivative of the matrix with respect to the parameters.
     *
     *  The first dimension of the array is the parameter dimension.
     */
    void evaluateModelDerivative(
        ndarray::Array<double,3,3> const & derivative,
        ndarray::Array<double const,1,1> const & param
    ) const;

    void writeInitialParameters(ndarray::Array<double,1,1> const & param) const;

    virtual ~Evaluator() {}

protected:

    BaseEvaluator(int data_size, int coefficient_size, int parameter_size) :
        _coefficient_size(coefficient_size),
        _parameter_size(parameter_size),
        _data_vector(ndarray::allocate(ndarray::makeVector(data_size)))
    {}

    BaseEvaluator(ndarray::Array<double,1,1> const & data, int coefficient_size, int parameter_size) :
        _coefficient_size(coefficient_size),
        _parameter_size(parameter_size),
        _data_vector(data)
    {}

    BaseEvaluator(Evaluator const & other) :
        _coefficient_size(other._coefficient_size), 
        _parameter_size(other._parameter_size), 
        _data_vector(other._data_vector)
    {}

    virtual void _evaluateModelMatrix(
        ndarray::Array<double,2,2> const & matrix,
        ndarray::Array<double const,1,1> const & param
    ) const = 0;

    virtual void _evaluateModelDerivative(
        ndarray::Array<double,3,3> const & derivative,
        ndarray::Array<double const,1,1> const & param
    ) const;

    virtual void _writeInitialParameters(ndarray::Array<double,1,1> const & param) const = 0;

    int const _coefficient_size;
    int const _parameter_size;
    ndarray::Array<double,1,1> _data_vector;

private:
    void operator=(BaseEvaluator const &) {}
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_BaseEvaluator
