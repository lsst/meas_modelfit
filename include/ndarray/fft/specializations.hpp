#ifndef NDARRAY_FFT_specializations_hpp_INCLUDED
#define NDARRAY_FFT_specializations_hpp_INCLUDED

/**
 *  @file ndarray/fft/specializations.hpp
 *  \brief Specializations of ArrayTraits, ExpressionTraits, and ArrayAccess for FourierArray.
 */

#include "ndarray.hpp"
#include "ndarray/fft_fwd.hpp"
#include "ndarray/fft/FourierTraits.hpp"

namespace ndarray {
namespace detail {

/**
 *  \internal \brief Common aspects of ArrayTraits for Array.
 *
 *  \ingroup FFTInternalGroup
 */
template <typename T, int N, int C>
struct FourierArrayTraitsBase {
    typedef boost::mpl::int_<N> ND;
    typedef boost::mpl::int_<C> RMC;

    template <typename T1>
    struct RebindT {
        typedef FourierArray<T1,N,C> Other;
    };

    template <int C1>
    struct RebindC {
        typedef FourierArray<T,N,C1> Other;
    };

    template <int N1, int C1>
    struct RebindN {
        typedef FourierArray<T,N1,C1> Other;
    };
};

/**
 *  \internal \brief ArrayTraits specialization for multidimensional FourierArray.
 *
 *  \ingroup FFTInternalGroup
 */
template <typename T, int N, int C>
struct ArrayTraits< FourierArray<T,N,C> > : public FourierArrayTraitsBase<T,N,C> {
    typedef typename detail::FourierTraits<T>::ElementK Element;
    typedef NestedIterator< FourierArray<T,N,C> > Iterator;
    typedef FourierArray<T,N-1,(N==C)?(N-1):C> Reference;
    typedef Reference Value;
};

/**
 *  \internal \brief ArrayTraits specialization for 1D contiguous FourierArray.
 *
 *  \ingroup FFTInternalGroup
 */
template <typename T>
struct ArrayTraits< FourierArray<T,1,1> > : public FourierArrayTraitsBase<T,1,1> {
    typedef typename detail::FourierTraits<T>::ElementK Element;
    typedef Element * Iterator;
    typedef Element & Reference;
    typedef Element Value;
};

/**
 *  \internal \brief ArrayTraits specialization for 1D noncontiguous FourierArray.
 *
 *  \ingroup FFTInternalGroup
 */
template <typename T>
struct ArrayTraits< FourierArray<T,1,0> > : public FourierArrayTraitsBase<T,1,0> {
    typedef typename detail::FourierTraits<T>::ElementK Element;
    typedef StridedIterator<Element> Iterator;
    typedef Element & Reference;
    typedef Element Value;
};

/**
 *  \internal \brief ExpressionTraits specialization for FourierArray.
 *
 *  \ingroup FFTInternalGroup
 */
template <typename T, int N, int C>
struct ExpressionTraits< FourierArray<T,N,C> > {
    typedef FourierArray<T,N,C> Self;
    typedef typename ArrayTraits<Self>::Element Element;
    typedef typename ArrayTraits<Self>::ND ND;
    typedef typename ArrayTraits<Self>::Iterator Iterator;
    typedef typename ArrayTraits<Self>::Reference Reference;
    typedef typename ArrayTraits<Self>::Value Value;
};

/**
 *  \internal
 *  \brief Specialization of ArrayAccess for FourierArray.
 *
 *  \ingroup FFTInternalGroup
 */
template <typename T, int N, int C>
struct ArrayAccess< FourierArray<T,N,C> > {
    typedef typename ArrayTraits< FourierArray<T,N,C> >::Element Element;
    typedef Array<Element,N,C> Base;

    static inline Element * & getArrayData(FourierArray<T,N,C> & array) {
        return ArrayAccess<Base>::getArrayData(array._base);
    }

    template <typename T1, int N1, int C1>
    static inline FourierArray<T,N,C> construct(
        Element * data, 
        FourierArray<T1,N1,C1> const & array
    ) {
        return FourierArray<T,N,C>(array._size,ArrayAccess<Base>::construct(data,array._base));
    }

    static inline void assign(FourierArray<T,N,C> & to, FourierArray<T,N,C> const & from) {
        ArrayAccess<Base>::assign(to._base,from._base);
        to._size = from._size;
    }
};

} // namespace ndarray::detail
} // namespace ndarray

#endif // !NDARRAY_FFT_specializations_hpp_INCLUDED
