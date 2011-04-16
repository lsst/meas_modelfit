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

#ifndef LSST_MEAS_MULTIFIT_SAMPLING_Table
#define LSST_MEAS_MULTIFIT_SAMPLING_Table

#include "lsst/ndarray.h"
#include "lsst/ndarray/eigen.h"
#include "lsst/meas/multifit/constants.h"
#include <boost/type_traits/is_convertible.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/noncopyable.hpp>

namespace lsst { namespace meas { namespace multifit { namespace sampling {

enum NestedSampleType {
    MACLAURIN_SERIES, // ln(f(0)), [d ln(f) / dx](0), [d^2 ln(f) / dx^2](0)
    FISHER_FULL,      // ln(f(x_0)), x_0, [d^2 ln(f) / dx^2](x_0)
    FISHER_LLT,       // like FISHER_FULL, but stores the lower Cholesky factor of the matrix
    COVARIANCE_FULL,  // like FISHER_FULL, but with the matrix inverted
    COVARIANCE_LLT,   // like COVARIANCE_FULL, but stores the lower Cholesky factor of the matrix
};

namespace detail {

template <typename T>
struct RecordBase : private boost::noncopyable {

    struct Nested : private boost::noncopyable {
        mutable T & scalar;
        ndarray::EigenView<T,1,1> vector;
        ndarray::EigenView<T,2,2> matrix;

        Nested(T & scalar_, ndarray::Array<T,1,1> const & vector_, ndarray::Array<T,2,2> const & matrix_) :
            scalar(scalar_), vector(vector_), matrix(matrix_) {}
    };

    mutable T & importance;
    mutable T & target;
    mutable T & weight;
    ndarray::ArrayRef<T,1,1> const parameters;
    mutable Nested nested;

protected:

    RecordBase(
        T & importance_, T & target_, T & weight_,
        lsst::ndarray::Array<T,1,1> const & parameters_,
        T & nestedScalar,
        lsst::ndarray::Array<T,1,1> const & nestedVector,
        lsst::ndarray::Array<T,2,2> const & nestedMatrix
    ) : importance(importance_), target(target_), weight(weight_), parameters(parameters_),
        nested(nestedScalar, nestedVector, nestedMatrix)
    {}

};

template <typename T>
struct TableBase : private boost::noncopyable {

    struct Nested : private boost::noncopyable {
        ndarray::ArrayRef<T,1,0> const scalar;
        ndarray::ArrayRef<T,2,1> const vector;
        ndarray::ArrayRef<T,3,2> const matrix;

        Nested(
            ndarray::Array<T,1,0> const & scalar_,
            ndarray::Array<T,2,1> const & vector_,
            ndarray::Array<T,3,2> const & matrix_
        ) : scalar(scalar_), vector(vector_), matrix(matrix_) {}
    };

    lsst::ndarray::ArrayRef<T,1,0> const importance;
    lsst::ndarray::ArrayRef<T,1,0> const target;
    lsst::ndarray::ArrayRef<T,1,0> const weight;
    lsst::ndarray::ArrayRef<T,2,1> const parameters;
    mutable Nested nested;

    int getTableSize() const { return importance.template getSize<0>(); }

    int getParameterSize() const { return parameters.template getSize<1>(); }

    int getNestedSize() const { return nested.vector.template getSize<1>(); }

    double getWeightSum() const { return _weightSum; }

    NestedSampleType getNestedSampleType() const { return _nestedSampleType; }
    
protected:

    TableBase(
        lsst::ndarray::Array<T,1,0> const & importance_,
        lsst::ndarray::Array<T,1,0> const & target_,
        lsst::ndarray::Array<T,1,0> const & weight_,
        lsst::ndarray::Array<T,2,1> const & parameters_,
        lsst::ndarray::Array<T,1,0> const & nestedScalar,
        lsst::ndarray::Array<T,2,1> const & nestedVector,
        lsst::ndarray::Array<T,3,2> const & nestedMatrix,
        double weightSum, NestedSampleType nestedSampleType
    ) : importance(importance_), target(target_), weight(weight_),
        parameters(parameters_),
        nested(nestedScalar, nestedVector, nestedMatrix),
        _weightSum(weightSum), _nestedSampleType(nestedSampleType)
    {}

    mutable T _weightSum;
    mutable NestedSampleType _nestedSampleType;
};

} // namespace detail

class Record : public detail::RecordBase<double> {
public:

    Record(Record const & other) :
        detail::RecordBase<double>(
            other.importance, other.target, other.weight, other.parameters,
            other.nested.scalar, other.nested.vector.getArray(), other.nested.matrix.getArray()
        )
    {}
    
    Record(detail::TableBase<double> const & table, int n) :
        detail::RecordBase<double>(
            table.importance[n], table.target[n], table.weight[n], table.parameters[n],
            table.nested.scalar[n], table.nested.vector[n], table.nested.matrix[n]
        )
    {}

};

class ConstRecord : public detail::RecordBase<double const> {
public:

    ConstRecord(ConstRecord const & other) :
        detail::RecordBase<double const>(
            other.importance, other.target, other.weight, other.parameters,
            other.nested.scalar, other.nested.vector.getArray(), other.nested.matrix.getArray()
        )
    {}

    ConstRecord(Record const & other) :
        detail::RecordBase<double const>(
            other.importance, other.target, other.weight, other.parameters,
            other.nested.scalar, other.nested.vector.getArray(), other.nested.matrix.getArray()
        )
    {}

    ConstRecord(detail::TableBase<double> const & table, int n) :
        detail::RecordBase<double const>(
            table.importance[n], table.target[n], table.weight[n], table.parameters[n],
            table.nested.scalar[n], table.nested.vector[n], table.nested.matrix[n]
        )
    {}

    ConstRecord(detail::TableBase<double const> const & table, int n) :
        detail::RecordBase<double const>(
            table.importance[n], table.target[n], table.weight[n], table.parameters[n],
            table.nested.scalar[n], table.nested.vector[n], table.nested.matrix[n]
        )
    {}
    
};

class Table : public detail::TableBase<double> {
public:

    Table(Table const & other) :
        detail::TableBase<double>(
            other.importance, other.target, other.weight, other.parameters,
            other.nested.scalar, other.nested.vector, other.nested.matrix,
            other._weightSum, other._nestedSampleType
        )
    {}

    static Table allocate(
        int tableSize, int parameterSize, 
        NestedSampleType nestedSampleType=MACLAURIN_SERIES, int nestedSize=0
    );

    void computeWeights() const;

    void clear() const;

    Record operator[](int n) const { return Record(*this, n); }

private:
    
    Table(
        lsst::ndarray::Array<double,1,0> const & importance_,
        lsst::ndarray::Array<double,1,0> const & target_,
        lsst::ndarray::Array<double,1,0> const & weight_,
        lsst::ndarray::Array<double,2,1> const & parameters_,
        lsst::ndarray::Array<double,1,0> const & nestedScalar,
        lsst::ndarray::Array<double,2,1> const & nestedVector,
        lsst::ndarray::Array<double,3,2> const & nestedMatrix,
        double weightSum, NestedSampleType nestedSampleType
    ) : 
        detail::TableBase<double>(
            importance_, target_, weight_, parameters_,
            nestedScalar, nestedVector, nestedMatrix,
            weightSum, nestedSampleType
        )
    {}

};

class ConstTable : public detail::TableBase<double const> {
public:

    ConstTable(ConstTable const & other) :
        detail::TableBase<double const>(
            other.importance, other.target, other.weight, other.parameters,
            other.nested.scalar, other.nested.vector, other.nested.matrix,
            other._weightSum, other._nestedSampleType
        )
    {}

    ConstTable(Table const & other) :
        detail::TableBase<double const>(
            other.importance, other.target, other.weight, other.parameters,
            other.nested.scalar, other.nested.vector, other.nested.matrix,
            other.getWeightSum(), other.getNestedSampleType()
        )
    {}

    ConstRecord operator[](int n) const { return ConstRecord(*this, n); }

};

}}}} // namespace lsst::meas::multifit::sampling

#endif // !LSST_MEAS_MULTIFIT_SAMPLING_Table
