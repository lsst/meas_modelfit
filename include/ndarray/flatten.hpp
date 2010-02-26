#ifndef NDARRAY_flatten_HPP_INCLUDED
#define NDARRAY_flatten_HPP_INCLUDED

#include "ndarray/Array.hpp"

namespace ndarray {

/**
 *  \ingroup MainGroup
 *  \brief Create a view into an array with trailing contiguous dimensions merged.
 *
 *  The first template parameter sets the dimension of the output array and must
 *  be specified directly.  Only row-major contiguous dimensions can be flattened.
 */
template <int Nf, typename T, int N, int C>
inline typename boost::enable_if_c< ((C+Nf-N)>=1), Array<T,Nf,(C+Nf-N)> >::type
flatten(Array<T,N,C> const & input) {
    BOOST_STATIC_ASSERT(C+Nf-N >= 1);
    Vector<int,N> oldShape = input.getShape();
    Vector<int,Nf> newShape = oldShape.template first<Nf>();
    for (int n=Nf; n<N; ++n)
        newShape[Nf-1] *= oldShape[n];
    Vector<int,Nf> newStrides = input.getStrides().template first<Nf>();
    newStrides[Nf-1] = 1;
    return ndarray::external(input.getData(), newShape, newStrides, input.getOwner());
}

} // namespace ndarray

#endif // !NDARRAY_flatten_HPP_INCLUDED
