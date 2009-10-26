#ifndef NDARRAY_DETAIL_ArrayTraits_hpp_INCLUDED
#define NDARRAY_DETAIL_ArrayTraits_hpp_INCLUDED

/** 
 *  @file ndarray/detail/ArrayTraits.hpp
 *
 *  \brief Traits for Array.
 */

#include <boost/type_traits/remove_const.hpp>
#include <boost/type_traits/add_const.hpp>

#include "ndarray_fwd.hpp"
#include "ndarray/detail/Core.hpp"

namespace ndarray {
namespace detail {

/**
 *  \internal \brief Traits for Array.
 *
 *  \ingroup InternalGroup
 */
template <typename ArrayT> struct ArrayTraits {};

/**
 *  \internal \brief Common aspects of ArrayTraits for Array.
 *
 *  \ingroup InternalGroup
 */
template <typename T, int N, int C>
struct ArrayTraitsBase {
    typedef T Element;
    typedef boost::mpl::int_<N> ND;
    typedef boost::mpl::int_<C> RMC;

    template <typename T1>
    struct RebindT {
        typedef Array<T1,N,C> Other;
    };

    template <int C1>
    struct RebindC {
        typedef Array<T,N,C1> Other;
    };

    template <int N1, int C1>
    struct RebindN {
        typedef Array<T,N1,C1> Other;
    };
};

/**
 *  \internal \brief ArrayTraits specialization for multidimensional Array.
 *
 *  \ingroup InternalGroup
 */
template <typename T, int N, int C>
struct ArrayTraits< Array<T,N,C> > : public ArrayTraitsBase<T,N,C> {
    typedef NestedIterator< Array<T,N,C> > Iterator;
    typedef Array<T,N-1,(N==C)?(N-1):C> Reference;
    typedef Reference Value;
};

/**
 *  \internal \brief ArrayTraits specialization for 1D contiguous Array.
 *
 *  \ingroup InternalGroup
 */
template <typename T>
struct ArrayTraits< Array<T,1,1> > : public ArrayTraitsBase<T,1,1> {
    typedef T * Iterator;
    typedef T & Reference;
    typedef T Value;
};

/**
 *  \internal \brief ArrayTraits specialization for 1D noncontiguous Array.
 *
 *  \ingroup InternalGroup
 */
template <typename T>
struct ArrayTraits< Array<T,1,0> > : public ArrayTraitsBase<T,1,0> {
    typedef StridedIterator<T> Iterator;
    typedef T & Reference;
    typedef T Value;
};

} // namespace ndarray::detail
} // namespace ndarray

#endif // !NDARRAY_DETAIL_ArrayTraits_hpp_INCLUDED
