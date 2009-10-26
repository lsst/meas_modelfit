#ifndef NDARRAY_FFT_FourierTransform_hpp_INCLUDED
#define NDARRAY_FFT_FourierTransform_hpp_INCLUDED

/** 
 *  @file ndarray/fft/FourierTransform.hpp
 *
 *  \brief Definitions for FourierTransform.
 */

#include <boost/noncopyable.hpp>

#include "ndarray.hpp"
#include "ndarray/fft/FourierArray.hpp"

namespace ndarray {

/**
 *  \ingroup FFTGroup
 *  \brief A wrapper for FFTW plans for fast Fourier transforms.
 *
 *  An instance of FourierTransform holds an FFTW "plan", providing repeated forward or
 *  inverse FFTs of predetermined arrays.
 *
 *  Multiplex plans can also be generated to perform an N-dimensional FFT on the nested arrays
 *  of an (N+1)-dimensional array.
 *
 *  Static member functions of FourierTransform are used to create instances, and optionally
 *  initialize the involved arrays.
 *
 */
template <typename T, int N>
class FourierTransform : private boost::noncopyable {
    BOOST_STATIC_ASSERT((!boost::is_const<T>::value));
public:
    
    typedef boost::shared_ptr<FourierTransform> Ptr;
    typedef boost::shared_ptr<FourierTransform const> ConstPtr;

    typedef Vector<int,N> Index; ///< Shape type for arrays.
    typedef Array<T,N,N> ArrayX; ///< Real-space array type.
    typedef FourierArray<T,N,N> ArrayK;  ///< Fourier-space array type.
    typedef Vector<int,N+1> MultiplexIndex; ///< Shape type for multiplexed arrays.
    typedef Array<T,N+1,N+1> MultiplexArrayX; ///< Real-space multiplexed array type.
    typedef FourierArray<T,N+1,N+1> MultiplexArrayK; ///< Fourier-space multiplexed array type.
    typedef typename ArrayX::Element ElementX; ///< Real-space array data type;
    typedef typename ArrayK::Element ElementK; ///< Fourier-space array data type;

    /**
     *  \brief Create a plan for forward-transforming a single N-dimensional array.
     *   
     *  Arrays will be initialized with new memory if empty.  If they are not empty,
     *  existing data may be overwritten when the plan is created.
     */
    static Ptr planForward(
        Index const & shape,  ///< Shape of the real-space array.
        ArrayX & x,           ///< Input real-space array.
        ArrayK & k            ///< Output Fourier-space array.
    );

    /**
     *  \brief Create a plan for inverse-transforming a single N-dimensional array.
     *   
     *  Arrays will be initialized with new memory if empty.  If they are not empty,
     *  existing data may be overwritten when the plan is created.
     */
    static Ptr planInverse(
        Index const & shape,  ///< Shape of the real-space array.
        ArrayK & k,           ///< Input Fourier-space array.
        ArrayX & x            ///< Output real-space array.
    );

    /**
     *  \brief Create a plan for forward-transforming a sequence of nested N-dimensional arrays.
     *   
     *  Arrays will be initialized with new memory if empty.  If they are not empty,
     *  existing data may be overwritten when the plan is created.
     */
    static Ptr planMultiplexForward(
        MultiplexIndex const & shape, ///< Shape of the real-space array. First dimension is multiplexed.
        MultiplexArrayX & x,          ///< Input real-space array.
        MultiplexArrayK & k           ///< Output Fourier-space array.
    );

    /**
     *  \brief Create a plan for inverse-transforming a sequence of nested N-dimensional arrays.
     *   
     *  Arrays will be initialized with new memory if empty.  If they are not empty,
     *  existing data may be overwritten when the plan is created.
     */
    static Ptr planMultiplexInverse(
        MultiplexIndex const & shape, ///< Shape of the real-space array. First dimension is multiplexed.
        MultiplexArrayK & k,          ///< Input Fourier-space array.
        MultiplexArrayX & x           ///< Output real-space array.
    );

    /// \brief Create a new real-space array with the given real-space shape.
    template <int M>
    static Array<T,M,M> initializeX(Vector<int,M> const & shape);

    /// \brief Create a new Fourier-space array with the given real-space shape.
    template <int M>
    static FourierArray<T,M,M> initializeK(Vector<int,M> const & shape);

    /** 
     *  \brief Initialize, as necessary, a pair of arrays with the given real-space shape.
     * 
     *  If either array is not empty, it must be consistent with the given shape.
     */
    template <int M>
    static void initialize(Vector<int,M> const & shape, Array<T,M,M> & x, FourierArray<T,M,M> & k);

    /// \brief Execute the FFTW plan.
    void execute();

    ~FourierTransform();

private:
    typedef boost::shared_ptr<ElementX> OwnerX;
    typedef boost::shared_ptr<ElementK> OwnerK;

    FourierTransform(void * plan, OwnerX const & x, OwnerK const & k)
        : _plan(plan), _x(x), _k(k) {}

    void * _plan; // 'void' so we don't have to include fftw3.h in the header file
    OwnerX _x;
    OwnerK _k;
};

} // namespace ndarray

#endif // !NDARRAY_FFT_FourierTransform_hpp_INCLUDED
