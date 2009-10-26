#ifndef NDARRAY_casts_hpp_INCLUDED
#define NDARRAY_casts_hpp_INCLUDED

/** 
 *  @file ndarray/casts.hpp
 *
 *  \brief Specialized casts for Array.
 */

#include "ndarray/Array.hpp"
#include <boost/type_traits/add_const.hpp>
#include <boost/type_traits/remove_const.hpp>

namespace ndarray {

/// \addtogroup MainGroup
/// @{

/**
 *  Convert an Array with a const data type to an array
 *  with a non-const data type.
 */
template <typename T, typename ArrayT>
typename detail::ArrayTraits<ArrayT>::template RebindT<T>::Other
const_array_cast(ArrayT const & array) {
    typedef typename detail::ArrayTraits<ArrayT>::template RebindT<T>::Other Output;
    return detail::ArrayAccess<Output>::construct(
        const_cast<typename Output::Element*>(array.getData()),
        array
    );
}

/**
 *  Convert an Array to a type with more guaranteed
 *  row-major-contiguous dimensions with no checking.
 */
template <int C, typename ArrayT>
typename detail::ArrayTraits<ArrayT>::template RebindC<C>::Other
static_array_cast(ArrayT const & array) {
    typedef typename detail::ArrayTraits<ArrayT>::template RebindC<C>::Other Output;
    return detail::ArrayAccess<Output>::construct(
        const_cast<typename Output::Element*>(array.getData()),
        array
    );
}

/**
 *  Convert an Array to a type with more guaranteed
 *  row-major-contiguous dimensions, if the strides
 *  of the array match the desired number of RMC
 *  dimensions.  If the cast fails, an empty Array
 *  is returned.
 */
template <int C, typename ArrayT>
typename detail::ArrayTraits<ArrayT>::template RebindC<C>::Other
dynamic_array_cast(ArrayT const & array) {
    typedef typename detail::ArrayTraits<ArrayT>::template RebindC<C>::Other Output;
    const static int N = Output::ND::value;
    Vector<int,N> shape = array.getShape();
    Vector<int,N> strides = array.getStrides();
    int n = 1;
    for (int i=1; i<=C; ++i) {
        if (strides[N-i] != n) return Output();
        n *= shape[N-i];
    }
    return static_array_cast<C>(array);
}

/// @}

}

#endif // !NDARRAY_casts_hpp_INCLUDED
