#ifndef NDARRAY_FFT_FFTWTraits_hpp_INCLUDED
#define NDARRAY_FFT_FFTWTraits_hpp_INCLUDED

/** 
 *  @file ndarray/fft/FFTWTraits.hpp
 *
 *  \brief Traits classes that wrap FFTW in a template-friendly interface.
 */

#include <complex>
#include <fftw3.h>
#include "ndarray/fft/FourierTraits.hpp"
    
#define NDARRAY_FFTW_TRAITS(TYPE,PREFIX)                                \
    template <> struct FFTWTraits<TYPE> {                               \
        BOOST_STATIC_ASSERT((!boost::is_const<TYPE>::value));           \
        typedef PREFIX ## _plan Plan;                                   \
        typedef FourierTraits<TYPE>::ElementX ElementX;                 \
        typedef FourierTraits<TYPE>::ElementK ElementK;                 \
        typedef boost::shared_ptr<ElementX> OwnerX;                     \
        typedef boost::shared_ptr<ElementK> OwnerK;                     \
        static inline Plan forward(int rank, const int *n, int howmany,	\
                                   ElementX *in, const int *inembed, int istride, int idist, \
                                   ElementK *out, const int *onembed, int ostride, int odist, \
                                   unsigned flags) {			\
            return PREFIX ## _plan_many_dft_r2c(rank,n,howmany,         \
                                                in,inembed,istride,idist, \
                                                reinterpret_cast<PREFIX ## _complex*>(out), \
                                                onembed,ostride,odist,  \
                                                flags);			\
        }                                                               \
        static inline Plan inverse(int rank, const int *n, int howmany, \
                                   ElementK *in, const int *inembed, int istride, int idist, \
                                   ElementX *out, const int *onembed, int ostride, int odist, \
                                   unsigned flags) {			\
            return PREFIX ## _plan_many_dft_c2r(rank,n,howmany,         \
                                                reinterpret_cast<PREFIX ## _complex*>(in), \
                                                inembed,istride,idist,  \
                                                out,onembed,ostride,odist, \
                                                flags);			\
        }                                                               \
        static inline void destroy(Plan p) { PREFIX ## _destroy_plan(p); } \
        static inline void execute(Plan p) { PREFIX ## _execute(p); }	\
        static inline OwnerX allocateX(int n) {                         \
            return OwnerX(                                              \
                reinterpret_cast<ElementX*>(                            \
                    PREFIX ## _malloc(sizeof(ElementX)*n)               \
                ),                                                      \
                PREFIX ## _free                                         \
            );                                                          \
        }                                                               \
        static inline OwnerK allocateK(int n) {                         \
            return OwnerK(                                              \
                reinterpret_cast<ElementK*>(                            \
                    PREFIX ## _malloc(sizeof(ElementK)*n)               \
                ),                                                      \
                PREFIX ## _free                                         \
            );                                                          \
        }                                                               \
    };                                                                  \
    template <> struct FFTWTraits< std::complex<TYPE> > {               \
        typedef PREFIX ## _plan Plan;                                   \
        typedef FourierTraits< std::complex<TYPE> >::ElementX ElementX; \
        typedef FourierTraits< std::complex<TYPE> >::ElementK ElementK; \
        typedef boost::shared_ptr<ElementX> OwnerX;                     \
        typedef boost::shared_ptr<ElementK> OwnerK;                     \
        static inline Plan forward(int rank, const int *n, int howmany,	\
                                   ElementX *in, const int *inembed, int istride, int idist, \
                                   ElementK *out, const int *onembed, int ostride, int odist, \
                                   unsigned flags) {			\
            return PREFIX ## _plan_many_dft(rank,n,howmany,             \
                                            reinterpret_cast<PREFIX ## _complex*>(in), \
                                            inembed,istride,idist,      \
                                            reinterpret_cast<PREFIX ## _complex*>(out), \
                                            onembed,ostride,odist,      \
                                            FFTW_FORWARD,flags);        \
        }                                                               \
        static inline Plan inverse(int rank, const int *n, int howmany, \
                                   ElementK *in, const int *inembed, int istride, int idist, \
                                   ElementX *out, const int *onembed, int ostride, int odist, \
                                   unsigned flags) {			\
            return PREFIX ## _plan_many_dft(rank,n,howmany,             \
                                            reinterpret_cast<PREFIX ## _complex*>(in), \
                                            inembed,istride,idist,      \
                                            reinterpret_cast<PREFIX ## _complex*>(out), \
                                            onembed,ostride,odist,      \
                                            FFTW_BACKWARD,flags);       \
        }                                                               \
        static inline void destroy(Plan p) { PREFIX ## _destroy_plan(p); } \
        static inline void execute(Plan p) { PREFIX ## _execute(p); }	\
        static inline OwnerX allocateX(int n) {                         \
            return OwnerX(                                              \
                reinterpret_cast<ElementX*>(                            \
                    PREFIX ## _malloc(sizeof(ElementX)*n)               \
                ),                                                      \
                PREFIX ## _free                                         \
            );                                                          \
        }                                                               \
        static inline OwnerK allocateK(int n) {                         \
            return OwnerK(                                              \
                reinterpret_cast<ElementK*>(                            \
                    PREFIX ## _malloc(sizeof(ElementK)*n)               \
                ),                                                      \
                PREFIX ## _free                                         \
            );                                                          \
        }                                                               \
    }

namespace ndarray {
namespace detail {

/**
 *  \internal \ingroup FFTInternalGroup
 *  \brief A traits class that maps C++ template types to FFTW types and wraps FFTW function calls.
 */
template <typename T> struct FFTWTraits { BOOST_STATIC_ASSERT(sizeof(T) < 0); };

/// \cond SPECIALIZATIONS

NDARRAY_FFTW_TRAITS(float,fftwf);
NDARRAY_FFTW_TRAITS(double,fftw);

#ifndef NDARRAY_NO_LONG_DOUBLE
NDARRAY_FFTW_TRAITS(long double, fftwl);
#endif

/// \endcond

} // namespace ndarray::detail
} // namespace ndarray

#undef NDARRAY_FFTW_TRAITS_REAL
#undef NDARRAY_FFTW_TRAITS_COMPLEX
#undef NDARRAY_FFTW_TRAITS

#endif // !NDARRAY_FFT_FFTWTraits_hpp_INCLUDED
