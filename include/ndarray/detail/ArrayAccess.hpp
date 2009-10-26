#ifndef NDARRAY_DETAIL_ArrayAccess_hpp_INCLUDED
#define NDARRAY_DETAIL_ArrayAccess_hpp_INCLUDED

/** 
 *  @file ndarray/detail/ArrayAccess.hpp
 *
 *  \brief Definitions for ArrayAccess
 */

#include "ndarray/Array.hpp"

namespace ndarray {
namespace detail {

/**
 *  \internal @class ArrayAccess
 *  \brief A friend class to Array that provides access to internals.
 *
 *  \ingroup InternalGroup
 *
 *  This intentionally partially breaks encapsulation for Array,
 *  allowing extension classes private access while hiding that
 *  access in a well-defined and clearly non-public interface.
 *
 *  ArrayAccess should be specialized for any type for which
 *  ArrayTraits is specialized; this will allow it to be used
 *  with the array cast operators, the shallow operator, and
 *  NestedIterator.
 *
 *  \todo Modify initializers and views for easier extension
 *  using this mechanism.
 */
template <typename ArrayT>
struct ArrayAccess {
    typedef typename ArrayTraits< ArrayT >::Element Element;

    static Element * & getArrayData(ArrayT & array);

    template <typename T1, int N1, int C1>
    static ArrayT construct(
        Element * data,
        typename ArrayT::template Rebind<T1,N1,C1>::Other const & array
    );

    static void assign(ArrayT & to, ArrayT const & from);
};

/**
 *  \internal
 *  \brief Specialization of ArrayAccess for Array.
 *
 *  \ingroup InternalGroup
 */
template <typename T, int N, int C>
struct ArrayAccess< Array<T,N,C> > {
    typedef typename ArrayTraits< Array<T,N,C> >::Element Element;

    static inline Element * & getArrayData(Array<T,N,C> & array) { return array._data; }

    template <typename T1, int N1, int C1>
    static inline Array<T,N,C> construct(Element * data, Array<T1,N1,C1> const & array) {
        return Array<T,N,C>(data,array._core);
    }

    static inline void assign(Array<T,N,C> & to, Array<T,N,C> const & from) {
        to._core = from._core;
        to._data = from._data;
    }
};

} // namespace ndarray::detail
} // namespace ndarray

#endif // !NDARRAY_DETAIL_ArrayAccess_hpp_INCLUDED
