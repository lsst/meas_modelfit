#ifndef NDARRAY_DETAIL_ExpressionTraits_hpp_INCLUDED
#define NDARRAY_DETAIL_ExpressionTraits_hpp_INCLUDED

/** 
 *  @file ndarray/detail/ExpressionTraits.hpp
 *
 *  \brief Traits for Expression.
 */

#include "ndarray/detail/ArrayTraits.hpp"
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>

namespace ndarray {
namespace detail {

/**
 *  \internal \brief Traits for Expression.
 *
 *  \ingroup InternalGroup
 */
template <typename ExpressionT> struct ExpressionTraits {
    typedef boost::mpl::true_ IsScalar;
};

/**
 *  \internal \brief ExpressionTraits specialization for Array.
 *
 *  \ingroup InternalGroup
 */
template <typename T, int N, int C>
struct ExpressionTraits< Array<T,N,C> > {
    typedef Array<T,N,C> Self;
    typedef typename ArrayTraits<Self>::Element Element;
    typedef typename ArrayTraits<Self>::ND ND;
    typedef typename ArrayTraits<Self>::Iterator Iterator;
    typedef typename ArrayTraits<Self>::Reference Reference;
    typedef typename ArrayTraits<Self>::Value Value;
    typedef boost::mpl::false_ IsScalar;
};

/**
 *  \internal \brief ExpressionTraits specialization for UnaryOpExpression.
 *
 *  \ingroup InternalGroup
 */
template <typename Operand, typename UnaryFunction, int N>
struct ExpressionTraits< UnaryOpExpression<Operand,UnaryFunction,N> > {
    typedef typename UnaryFunction::result_type Element;
    typedef typename ExpressionTraits<Operand>::ND ND;
    typedef UnaryOpIterator<Operand,UnaryFunction> Iterator;
    typedef UnaryOpExpression<
        typename ExpressionTraits<Operand>::Reference,UnaryFunction,N-1
        > Value;
    typedef Value Reference;
    typedef boost::mpl::false_ IsScalar;
};

/**
 *  \internal \brief ExpressionTraits specialization for 1D UnaryOpExpression.
 *
 *  \ingroup InternalGroup
 */
template <typename Operand, typename UnaryFunction>
struct ExpressionTraits< UnaryOpExpression<Operand,UnaryFunction,1> > {
    typedef typename UnaryFunction::result_type Element;
    typedef typename ExpressionTraits<Operand>::ND ND;
    typedef UnaryOpIterator<Operand,UnaryFunction> Iterator;
    typedef Element Reference;
    typedef Element Value;
    typedef boost::mpl::false_ IsScalar;
};

/**
 *  \internal \brief ExpressionTraits specialization for BinaryOpExpression.
 *
 *  \ingroup InternalGroup
 */
template <typename Operand1, typename Operand2, typename BinaryFunction, int N>
struct ExpressionTraits< BinaryOpExpression<Operand1,Operand2,BinaryFunction,N> > {
    typedef typename BinaryFunction::result_type Element;
    typedef typename ExpressionTraits<Operand1>::ND ND;
    typedef BinaryOpIterator<Operand1,Operand2,BinaryFunction> Iterator;
    typedef BinaryOpExpression<
        typename ExpressionTraits<Operand1>::Reference,
        typename ExpressionTraits<Operand2>::Reference,
        BinaryFunction, N-1 > Reference;
    typedef Reference Value;
    typedef boost::mpl::false_ IsScalar;
    BOOST_STATIC_ASSERT((ND::value == ExpressionTraits<Operand2>::ND::value));
};

/**
 *  \internal \brief ExpressionTraits specialization for 1D BinaryOpExpression.
 *
 *  \ingroup InternalGroup
 */
template <typename Operand1, typename Operand2, typename BinaryFunction>
struct ExpressionTraits< BinaryOpExpression<Operand1,Operand2,BinaryFunction,1> > {
    typedef typename BinaryFunction::result_type Element;
    typedef typename ExpressionTraits<Operand1>::ND ND;
    typedef BinaryOpIterator<Operand1,Operand2,BinaryFunction> Iterator;
    typedef Element Value;
    typedef Element Reference;
    typedef boost::mpl::false_ IsScalar;
    BOOST_STATIC_ASSERT((ND::value == ExpressionTraits<Operand2>::ND::value));
};

} // namespace ndarray::detail
} // namespace ndarray

#endif // !NDARRAY_DETAIL_ExpressionTraits_hpp_INCLUDED
