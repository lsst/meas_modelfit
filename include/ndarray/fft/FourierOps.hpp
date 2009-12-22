#ifndef NDARRAY_FFT_FourierOps_hpp_INCLUDED
#define NDARRAY_FFT_FourierOps_hpp_INCLUDED

/** 
 *  @file ndarray/fft/FourierOps.hpp
 *
 *  \brief Definitions for FourierOps.
 */

#include <boost/noncopyable.hpp>

#include "ndarray.hpp"
#include "ndarray/fft/FourierArray.hpp"

namespace ndarray {
/// \cond INTERNAL
namespace detail {

/**
 *  \internal \ingroup FFTInternalGroup
 *  \brief Implementations for the shift() and differentiate() functions.
 */
template <typename T, int N>
struct FourierOps {
    typedef typename FourierTraits<T>::ElementK ElementK;
    
    template <int C>
    static void shift(T const * offset, ElementK const & factor, FourierArray<T,N,C> const & array) {
        typename FourierArray<T,N,C>::Iterator iter = array.begin();
        T u = -2.0 * M_PI * (*offset) / array.size();
        int kMid = (array.size() + 1) / 2;
        for (int k = 0; k < kMid; ++k, ++iter) {
            FourierOps<T,N-1>::shift(offset+1, factor * std::polar(static_cast<T>(1), u * k), *iter);
        }
        if (array.size() % 2 == 0) {
            FourierOps<T,N-1>::shift(offset+1, factor * std::cos(u * kMid), *iter);
            ++iter;
            ++kMid;
        }
        for (int k_n = kMid - array.size(); k_n < 0; ++k_n, ++iter) {
            FourierOps<T,N-1>::shift(offset+1, factor * std::polar(static_cast<T>(1), u * k_n), *iter);
        }
    }

    template <int C>
    static void differentiate(int m, FourierArray<T,N,C> const & array) {
        typename FourierArray<T,N,C>::Iterator iter = array.begin();
        int kMid = (array.size() + 1) / 2;
        T u = 2.0 * M_PI / array.size();
        for (int k = 0; k < kMid; ++k, ++iter) {
            if (m == N) (*iter) *= ElementK(static_cast<T>(0),u * k);
            FourierOps<T,N-1>::differentiate(m,*iter);
        }
        if (array.size() % 2 == 0) {
            (*iter) = static_cast<T>(0);
            ++iter;
            ++kMid;
        }
        for (int k_n = kMid - array.size(); k_n < 0; ++k_n, ++iter) {
            if (m == N) (*iter) *= ElementK(static_cast<T>(0),u * k_n);
            FourierOps<T,N-1>::differentiate(m,*iter);
        }
    }

};

/**
 *  \internal \ingroup FFTInternalGroup
 *  \brief Implementations for the shift() and differentiate() functions (1d specialization).
 */
template <typename T>
struct FourierOps<T,1> {
    typedef typename FourierTraits<T>::ElementK ElementK;
    
    template <int C>
    static void shift(T const * offset, ElementK const & factor, FourierArray<T,1,C> const & array) {
        typename FourierArray<T,1,C>::Iterator iter = array.begin();
        T u = -2.0 * M_PI * (*offset) / array.getRealSize();
        int kMid = (array.getRealSize() + 1) / 2;
        for (int k = 0; k < kMid; ++k, ++iter) {
            (*iter) *= factor * std::polar(1.0, u * k);
        }
        if (array.getRealSize() % 2 == 0) {
            (*iter) *= factor * std::cos(u * kMid);
            ++iter;
        }
    }

    template <int C>
    static void differentiate(int m, FourierArray<T,1,C> const & array) {
        typename FourierArray<T,1,C>::Iterator iter = array.begin();
        int kMid = (array.getRealSize() + 1) / 2;
        if (m == 1) {
            T u = 2.0 * M_PI / array.getRealSize();
            for (int k = 0; k < kMid; ++k, ++iter) {
                (*iter) *= ElementK(static_cast<T>(0),u * k);
            }
        }
        if (array.getRealSize() % 2 == 0) {
            array[kMid] = static_cast<T>(0);
        }            
    }

};

} // namespace ndarray::detail
/// \endcond

/**
 *  \brief Perform a Fourier-space translation transform.
 *
 *  \ingroup FFTGroup
 */
template <typename T, int N, int C>
void shift(Vector<T,N> const & offset, FourierArray<T,N,C> const & array) {
    detail::FourierOps<T,N>::shift(
        offset.begin(),
        static_cast<typename detail::FourierTraits<T>::ElementK>(1),
        array
    );
}

/**
 *  \brief Numerically differentiate the array in Fourier-space in the given dimension.
 *
 *  \ingroup FFTGroup
 */
template <typename T, int N, int C>
void differentiate(int n, FourierArray<T,N,C> const & array) {
    detail::FourierOps<T,N>::differentiate(N-n,array);
}

} // namespace ndarray

#endif // !NDARRAY_FFT_FourierOps_hpp_INCLUDED
