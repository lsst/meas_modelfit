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

#define NDARRAY_FUNCTION_TAG(NAME,BASE,OP)                              \
    /** \internal @class NAME  */                                       \
    /** \brief An AdaptableFunctionTag for the OP operation. */ \
    /** \ingroup InternalGroup */                                       \
    struct NAME : public AdaptableFunctionTag<NAME> {                   \
        template <typename A, typename B>                               \
        struct ScalarFunction : public BASE<A,B> {                      \
            typename BASE<A,B>::result_type operator()(                 \
                typename BASE<A,B>::ParamA a,                           \
                typename BASE<A,B>::ParamB b                            \
            ) const {                                                   \
                return a OP b;                                          \
            }                                                           \
        };                                                              \
    };
#ifndef DOXYGEN
#define NDARRAY_BINARY_OP(TAG,OP)                                       \
    /** \brief A lazy 'array OP scalar' operator. */                    \
    template <typename Operand, typename Scalar>                        \
    typename boost::enable_if<                                          \
        typename detail::ExpressionTraits<Scalar>::IsScalar,            \
        detail::UnaryOpExpression< Operand, typename TAG::template ExprScalar<Operand,Scalar>::Bound > \
    >::type                                                              \
    operator OP(Expression<Operand> const & operand, Scalar const & scalar) { \
        return vectorize(TAG::template ExprScalar<Operand,Scalar>::bind(scalar),operand); \
    }                                                                   \
    /** \brief A lazy 'scalar OP array' operator. */                    \
    template <typename Operand, typename Scalar>                        \
    typename boost::enable_if<                                          \
        typename detail::ExpressionTraits<Scalar>::IsScalar,            \
        detail::UnaryOpExpression< Operand, typename TAG::template ScalarExpr<Operand,Scalar>::Bound > \
    >::type                                                              \
    operator OP(Scalar const & scalar, Expression<Operand> const & operand) { \
        return vectorize(TAG::template ScalarExpr<Operand,Scalar>::bind(scalar),operand); \
    }                                                                   \
    /** \brief A lazy 'array OP array' operator. */                    \
    template <typename Operand1, typename Operand2>                     \
    detail::BinaryOpExpression< Operand1, Operand2,                     \
                                typename TAG::template ExprExpr<Operand1,Operand2>::BinaryFunction > \
    operator OP(Expression<Operand1> const & operand1, Expression<Operand2> const & operand2) { \
        return vectorize(                                               \
            typename TAG::template ExprExpr<Operand1,Operand2>::BinaryFunction(), \
            operand1,                                                   \
            operand2                                                   \
        );                                                              \
    }
#define NDARRAY_UNARY_OP(FUNCTOR,OP)                                    \
    /** \brief A lazy unary OP operator for arrays. */                  \
    template <typename Operand>                                         \
    detail::UnaryOpExpression< Operand, FUNCTOR<typename detail::ExpressionTraits<Operand>::Element> > \
    operator OP(Expression<Operand> const & operand) {                  \
        return vectorize(FUNCTOR<typename detail::ExpressionTraits<Operand>::Element>(),operand); \
    }
#else // DOXYGEN
#define NDARRAY_BINARY_OP(TAG,OP)                                       \
    /** \brief A lazy 'array OP scalar' operator. */                    \
    template <typename Operand, typename Scalar>                        \
    <unspecified-expression-type>                                      \
    operator OP(Expression<Operand> const & operand, Scalar const & scalar); \
    /** \brief A lazy 'scalar OP array' operator. */                    \
    template <typename Operand, typename Scalar>                        \
    <unspecified-expression-type>                                       \
    operator OP(Scalar const & scalar, Expression<Operand> const & operand); \
    /** \brief A lazy 'array OP array' operator. */                    \
    template <typename Operand1, typename Operand2>                     \
    <unspecified-expression-type>                                       \
    operator OP(Expression<Operand1> const & operand1, Expression<Operand2> const & operand2);
#define NDARRAY_UNARY_OP(FUNCTOR,OP)                                    \
    /** \brief A lazy unary OP operator for arrays. */                  \
    template <typename Operand>                                         \
    <unspecified-expression-type>                                       \
    operator OP(Expression<Operand> const & operand);
#endif

namespace ndarray {

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

NDARRAY_FUNCTION_TAG(PlusTag,PromotingBinaryFunction,+)
NDARRAY_FUNCTION_TAG(MinusTag,PromotingBinaryFunction,-)
NDARRAY_FUNCTION_TAG(MultipliesTag,PromotingBinaryFunction,*)
NDARRAY_FUNCTION_TAG(DividesTag,PromotingBinaryFunction,/)
NDARRAY_FUNCTION_TAG(ModulusTag,PromotingBinaryFunction,%)
NDARRAY_FUNCTION_TAG(BitwiseXorTag,PromotingBinaryFunction,^)
NDARRAY_FUNCTION_TAG(BitwiseOrTag,PromotingBinaryFunction,|)
NDARRAY_FUNCTION_TAG(BitwiseAndTag,PromotingBinaryFunction,|)
NDARRAY_FUNCTION_TAG(BitwiseLeftShiftTag,PromotingBinaryFunction,<<)
NDARRAY_FUNCTION_TAG(BitwiseRightShiftTag,PromotingBinaryFunction,>>)

NDARRAY_FUNCTION_TAG(EqualToTag,BinaryPredicate,==)
NDARRAY_FUNCTION_TAG(NotEqualToTag,BinaryPredicate,!=)
NDARRAY_FUNCTION_TAG(LessTag,BinaryPredicate,<)
NDARRAY_FUNCTION_TAG(GreaterTag,BinaryPredicate,>)
NDARRAY_FUNCTION_TAG(LessEqualTag,BinaryPredicate,<=)
NDARRAY_FUNCTION_TAG(GreaterEqualTag,BinaryPredicate,>=)
NDARRAY_FUNCTION_TAG(LogicalAnd,BinaryPredicate,&&)
NDARRAY_FUNCTION_TAG(LogicalOr,BinaryPredicate,||)

} // namespace ndarray::detail

/// \addtogroup OpGroup
/// @{
NDARRAY_BINARY_OP(detail::PlusTag,+)
NDARRAY_BINARY_OP(detail::MinusTag,-)
NDARRAY_BINARY_OP(detail::MultipliesTag,*)
NDARRAY_BINARY_OP(detail::DividesTag,/)
NDARRAY_BINARY_OP(detail::ModulusTag,%)
NDARRAY_BINARY_OP(detail::BitwiseXorTag,^)
NDARRAY_BINARY_OP(detail::BitwiseOrTag,|)
NDARRAY_BINARY_OP(detail::BitwiseAndTag,|)
NDARRAY_BINARY_OP(detail::BitwiseLeftShiftTag,<<)
NDARRAY_BINARY_OP(detail::BitwiseRightShiftTag,>>)

NDARRAY_BINARY_OP(detail::EqualToTag,==)
NDARRAY_BINARY_OP(detail::NotEqualToTag,!=)
NDARRAY_BINARY_OP(detail::LessTag,<)
NDARRAY_BINARY_OP(detail::GreaterTag,>)
NDARRAY_BINARY_OP(detail::LessEqualTag,<=)
NDARRAY_BINARY_OP(detail::GreaterEqualTag,>=)
NDARRAY_BINARY_OP(detail::LogicalAnd,&&)
NDARRAY_BINARY_OP(detail::LogicalOr,||)

NDARRAY_UNARY_OP(std::negate,-)
NDARRAY_UNARY_OP(std::logical_not,!)
NDARRAY_UNARY_OP(detail::BitwiseNot,~)
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

#undef NDARRAY_FUNCTION_TAG
#undef NDARRAY_BINARY_OP
#undef NDARRAY_UNARY_OP

#endif // !NDARRAY_operators_hpp_INCLUDED
