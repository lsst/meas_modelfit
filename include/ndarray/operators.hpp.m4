include(`operators.macros.m4')dnl
changecom(`###')dnl
#ifndef NDARRAY_operators_hpp_INCLUDED
#define NDARRAY_operators_hpp_INCLUDED

/** 
 *  \file ndarray/operators.hpp \brief Arithmetic and logical operators for Array.
 */

#include "ndarray/Array.hpp"
#include <boost/call_traits.hpp>
#include <boost/functional.hpp>

#include "ndarray/detail/UnaryOp.hpp"
#include "ndarray/detail/BinaryOp.hpp"

namespace ndarray {

/// \cond INTERNAL
namespace detail {

/** 
 *  \internal @class PromotingBinaryFunction
 *  \brief A typedef-providing base class for binary functors with numeric type promotion.
 *
 *  \ingroup InternalGroup
 */
template <typename A, typename B>
struct PromotingBinaryFunction {
    typedef A ElementA;
    typedef B ElementB;
    typedef A first_argument_type;
    typedef B second_argument_type;
    typedef typename boost::call_traits<A>::param_type ParamA;
    typedef typename boost::call_traits<B>::param_type ParamB;
    typedef typename Promote<A,B>::Type result_type;
};

/** 
 *  \internal @class BinaryPredicate
 *  \brief A typedef-providing base class for binary predicates.
 *
 *  \ingroup InternalGroup
 */
template <typename A, typename B>
struct BinaryPredicate {
    typedef A ElementA;
    typedef B ElementB;
    typedef A first_argument_type;
    typedef B second_argument_type;
    typedef typename boost::call_traits<A>::param_type ParamA;
    typedef typename boost::call_traits<B>::param_type ParamB;
    typedef bool result_type;
};

/** 
 *  \internal @class AdaptableFunctionTag
 *  \brief A CRTP base class for non-template classes that contain a templated functor.
 *
 *  \ingroup InternalGroup
 */
template <typename Derived>
struct AdaptableFunctionTag {

    template <typename OperandB, typename A>
    struct ScalarExpr {
        typedef typename Derived::template ScalarFunction<
            A, typename detail::ExpressionTraits<OperandB>::Element
            > BinaryFunction;
        typedef boost::binder1st<BinaryFunction> Bound;
        static Bound bind(A const & scalar) {
            return Bound(BinaryFunction(),scalar);
        }
    };

    template <typename OperandA, typename B>
    struct ExprScalar {
        typedef typename Derived::template ScalarFunction<
            typename detail::ExpressionTraits<OperandA>::Element, B
            > BinaryFunction;
        typedef boost::binder2nd<BinaryFunction> Bound;
        static Bound bind(B const & scalar) {
            return Bound(BinaryFunction(),scalar);
        }
    };

    template <typename OperandA, typename OperandB>
    struct ExprExpr {
        typedef typename Derived::template ScalarFunction<
            typename detail::ExpressionTraits<OperandA>::Element,
            typename detail::ExpressionTraits<OperandB>::Element
            > BinaryFunction;
    };

};

/**
 *  \internal @class BitwiseNot
 *  \ingroup InternalGroup
 *  \brief An STL Unary Function class for bitwise NOT (unary ~).
 */
template <typename T>
struct BitwiseNot {
    typedef T argument_type;
    typedef T result_type;

    result_type operator()(argument_type arg) const { return ~arg; }
};

FUNCTION_TAG(PlusTag,PromotingBinaryFunction,+)
FUNCTION_TAG(MinusTag,PromotingBinaryFunction,-)
FUNCTION_TAG(MultipliesTag,PromotingBinaryFunction,*)
FUNCTION_TAG(DividesTag,PromotingBinaryFunction,/)
FUNCTION_TAG(ModulusTag,PromotingBinaryFunction,%)
FUNCTION_TAG(BitwiseXorTag,PromotingBinaryFunction,^)
FUNCTION_TAG(BitwiseOrTag,PromotingBinaryFunction,|)
FUNCTION_TAG(BitwiseAndTag,PromotingBinaryFunction,|)
FUNCTION_TAG(BitwiseLeftShiftTag,PromotingBinaryFunction,<<)
FUNCTION_TAG(BitwiseRightShiftTag,PromotingBinaryFunction,>>)

FUNCTION_TAG(EqualToTag,BinaryPredicate,==)
FUNCTION_TAG(NotEqualToTag,BinaryPredicate,!=)
FUNCTION_TAG(LessTag,BinaryPredicate,<)
FUNCTION_TAG(GreaterTag,BinaryPredicate,>)
FUNCTION_TAG(LessEqualTag,BinaryPredicate,<=)
FUNCTION_TAG(GreaterEqualTag,BinaryPredicate,>=)
FUNCTION_TAG(LogicalAnd,BinaryPredicate,&&)
FUNCTION_TAG(LogicalOr,BinaryPredicate,||)

} // namespace ndarray::detail
/// \endcond

/// \addtogroup OpGroup
/// @{
BINARY_OP(detail::PlusTag,+)
BINARY_OP(detail::MinusTag,-)
BINARY_OP(detail::MultipliesTag,*)
BINARY_OP(detail::DividesTag,/)
BINARY_OP(detail::ModulusTag,%)
BINARY_OP(detail::BitwiseXorTag,^)
BINARY_OP(detail::BitwiseOrTag,|)
BINARY_OP(detail::BitwiseAndTag,|)
BINARY_OP(detail::BitwiseLeftShiftTag,<<)
BINARY_OP(detail::BitwiseRightShiftTag,>>)

BINARY_OP(detail::EqualToTag,==)
BINARY_OP(detail::NotEqualToTag,!=)
BINARY_OP(detail::LessTag,<)
BINARY_OP(detail::GreaterTag,>)
BINARY_OP(detail::LessEqualTag,<=)
BINARY_OP(detail::GreaterEqualTag,>=)
BINARY_OP(detail::LogicalAnd,&&)
BINARY_OP(detail::LogicalOr,||)

UNARY_OP(std::negate,-)
UNARY_OP(std::logical_not,!)
UNARY_OP(detail::BitwiseNot,~)
/// @}

template <typename Scalar>
inline typename boost::enable_if<typename detail::ExpressionTraits<Scalar>::IsScalar, bool>::type
any(Scalar const & scalar) {
    return bool(scalar);
}

/**
 *  \brief Return true if any of the elements of the given expression are true.
 *
 *  \ingroup MainGroup
 */
template <typename Derived>
inline bool
any(Expression<Derived> const & expr) {
    typename Derived::Iterator const i_end = expr.end();
    for (typename Derived::Iterator i = expr.begin(); i != i_end; ++i) {
        if (any(*i)) return true;
    }
    return false;
}

template <typename Scalar>
inline typename boost::enable_if<typename detail::ExpressionTraits<Scalar>::IsScalar, bool>::type
all(Scalar const & scalar) {
    return bool(scalar);
}

/**
 *  \brief Return true if all of the elements of the given expression are true.
 *
 *  \ingroup MainGroup
 */
template <typename Derived>
inline bool
all(Expression<Derived> const & expr) {
    typename Derived::Iterator const i_end = expr.end();
    for (typename Derived::Iterator i = expr.begin(); i != i_end; ++i) {
        if (!all(*i)) return false;
    }
    return true;
}

template <typename Scalar>
inline typename boost::enable_if<typename detail::ExpressionTraits<Scalar>::IsScalar, Scalar>::type
sum(Scalar const & scalar) { return scalar; }


/**
 *  \brief Return the sum of all elements of the given expression.
 *
 *  \ingroup MainGroup
 */
template <typename Derived>
inline typename Derived::Element
sum(Expression<Derived> const & expr) {
    typename Derived::Iterator const i_end = expr.end();
    typename Derived::Element total = static_cast<typename Derived::Element>(0);
    for (typename Derived::Iterator i = expr.begin(); i != i_end; ++i) {
        total += sum(*i);
    }
    return total;
}


} // namespace ndarray

#endif // !NDARRAY_operators_hpp_INCLUDED
