#ifndef NDARRAY_FFT_FourierArray_hpp_INCLUDED
#define NDARRAY_FFT_FourierArray_hpp_INCLUDED

/** 
 *  @file ndarray/fft/FourierArray.hpp
 *
 *  \brief Declarations for FourierArray.
 */

#include "ndarray.hpp"
#include "ndarray/fft/specializations.hpp"

#ifndef DOXYGEN
#define NDARRAY_ASSIGN(OP)                                              \
    template <typename OtherT>                                          \
    FourierArray const &                                                \
    operator OP(Expression<OtherT> const & expr) const {                \
        _base OP expr;                                                  \
        return *this;                                                   \
    }                                                                   \
    template <typename ScalarT>                                         \
    typename boost::enable_if<boost::is_convertible<ScalarT,Element>, FourierArray const &>::type \
    operator OP(ScalarT const & scalar) const {                         \
        _base OP scalar;                                                \
        return *this;                                                   \
    }
#else // DOXYGEN
#define NDARRAY_ASSIGN(OP)                                              \
    template <typename OtherT>                                          \
    FourierArray const &                                                \
    operator OP(Expression<OtherT> const & expr) const;                 \
    template <typename ScalarT>                                         \
    FourierArray const &                                                \
    operator OP(ScalarT const & scalar) const;
#endif

namespace ndarray {

/**
 *  \ingroup FFTGroup
 *  \brief A specialized Array-like class for Fourier-space arrays.
 *
 *  FourierArray provides almost all of the functionality of Array, but also knows
 *  the size of its real-space counterpart (which may be different for real-data
 *  transforms).
 *
 *  FourierArray objects can be most easily constructed using FourierTransform::initializeK().
 */
template <typename T, int N, int C>
class FourierArray : public detail::ArrayImpl< FourierArray<T,N,C> > {
public:
    /// \brief Data type of array elements.
    typedef typename detail::ArrayTraits<FourierArray>::Element Element;
    /// \brief Number of dimensions (boost::mpl::int_).
    typedef typename detail::ArrayTraits<FourierArray>::ND ND;
    /// \brief Number of guaranteed row-major contiguous dimensions, counted from the end (boost::mpl::int_).
    typedef typename detail::ArrayTraits<FourierArray>::RMC RMC;
    /// \brief Nested array or element iterator.
    typedef typename detail::ArrayTraits<FourierArray>::Iterator Iterator; 
    /// \brief Nested array or element reference.
    typedef typename detail::ArrayTraits<FourierArray>::Reference Reference;
    /// \brief Nested array or element value type.
    typedef typename detail::ArrayTraits<FourierArray>::Value Value;
    /// \brief Shared pointer that manages memory.
    typedef boost::shared_ptr<Element> Owner;
    /// \brief Vector type for N-dimensional indices.
    typedef Vector<int,N> Index;
    /// \brief Underlying array type.
    typedef Array<Element,N,C> Base;

    /** 
     *  \brief Default constructor. 
     *
     *  Creates an empty array with zero dimensions and null memory.
     */
    FourierArray() : _size(0), _base() {}

    /** 
     *  \brief Converting copy constructor. 
     *
     *  Implicit conversion is allowed for non-const -> const and for 
     *  more guaranteed RMC -> less guaranteed RMC (see \ref overview).
     */
    template <typename T1, int C1>
    FourierArray(
        FourierArray<T1,N,C1> const & other
#ifndef DOXYGEN
        , typename boost::enable_if_c<((C1>=C) && boost::is_convertible<T1*,T*>::value),void*>::type=0
#endif
    ) : _size(other._size), _base(other._base)  {}

    /// \brief Return a raw pointer to the beginning of the array.
    Element * getData() const { return this->_base.getData(); }

    /// \brief Return the shared_ptr that manages the lifetime of the array data.
    Owner getOwner() const { return this->_base.getOwner(); }

    /**
     *  \brief Return the size of a specific dimension.
     *
     *  The size of the last dimension is the size of the underlying complex array,
     *  not the size of the real-space array.
     */
    template <int P> int getSize() const { return this->_base.template getSize<P>(); }

    /// \brief Return the stride in a specific dimension.
    template <int P> int getStride() const { return this->_base.template getStride<P>(); }

    /**
     *  \brief Return a Vector of the sizes of all dimensions.
     *
     *  The size of the last dimension is the size of the underlying complex array,
     *  not the size of the real-space array.
     */
    Index getShape() const { return this->_base.getShape(); }

    /// \brief Return a Vector of the strides of all dimensions.
    Index getStrides() const { return this->_base.getStrides(); }

    /// \briref Return the total number of elements in the array.
    int getNumElements() const { return this->_base.getNumElements(); }

    /// \internal \brief A template metafunction class to determine the result of a view indexing operation.
    template <
        typename ViewT, 
        typename SeqT = typename detail::ViewNormalizer<N,typename ViewT::Sequence>::Output
        >
    struct ResultOf {
        typedef T Element;
        typedef typename detail::ViewTraits<N,C,SeqT>::ND ND;
        typedef typename detail::ViewTraits<N,C,SeqT>::RMC RMC;
        typedef FourierArray<Element,ND::value,RMC::value> Type;

        BOOST_STATIC_ASSERT((boost::is_same<typename boost::fusion::result_of::back<SeqT>::type,
                             detail::FullDim&>::value));
    };

    /// \brief Return a general view into this array (see \ref tutorial).
    template <typename ViewT>
    typename ResultOf<ViewT>::Type
    operator[](ViewT const & def) const {
        return typename ResultOf<ViewT>::Type(_size,_base[def]);
    }

    /// \brief Return a single subarray (for ND > 1) or element for (ND==1).
    using detail::ArrayImpl<FourierArray>::operator[];

    /// \brief Deep assignment operator.
    FourierArray const & operator=(FourierArray const & other) const {
        if (&other != this) {
            NDARRAY_ASSERT(other.getShape() == this->getShape());
            NDARRAY_ASSERT(other._size == this->_size);
            std::copy(other.begin(),other.end(),this->begin());
        }
        return *this;
    }

    /// \brief Return the underlying complex Array object.
    Base const & getBase() const { return _base; }

    /// \brief Return the size of the last dimension of the corresponding real-space array.
    int getRealSize() const { return _size; }

    /// \brief Return the shape of the corresponding real-space array.
    Index getRealShape() const {
        Index r = this->_base.getShape();
        r[N-1] = getRealSize();
        return r; 
    }

    /**
     *  \brief Construct a FourierArray from the real-space size in the last dimension and a Base array.
     *
     *  Note that this also allows construction via the allocate() and external() initializers, if
     *  their return values are passed as the second argument of the FourierArray constructor.
     */
    FourierArray(int size, Base const & base) :
        _size(size), _base(base) {
        NDARRAY_ASSERT((detail::FourierTraits<T>::computeLastDimensionSize(size) 
                        == _base.template getSize<N-1>()));
    }

    NDARRAY_ASSIGN(=)
    NDARRAY_ASSIGN(+=)
    NDARRAY_ASSIGN(-=)
    NDARRAY_ASSIGN(*=)
    NDARRAY_ASSIGN(/=)

private:
    template <typename T1, int N1, int C1> friend class FourierArray;
    template <typename ArrayT1> friend class detail::ArrayAccess;

    int _size;
    Base _base;
};

} // namespace ndarray

#undef NDARRAY_ASSIGN

#endif // !NDARRAY_FFT_FourierArray_hpp_INCLUDED
