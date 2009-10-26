#ifndef NDARRAY_Array_hpp_INCLUDED
#define NDARRAY_Array_hpp_INCLUDED

/** 
 *  @file ndarray/Array.hpp
 *
 *  \brief Definitions for Array.
 */

#include "ndarray/Vector.hpp"
#include "ndarray/detail/Core.hpp"
#include "ndarray/detail/ArrayTraits.hpp"
#include "ndarray/detail/ArrayImpl.hpp"
#include "ndarray/detail/ArrayAccess.hpp"
#include "ndarray/views.hpp"

#ifndef DOXYGEN
#define NDARRAY_GENERAL_ASSIGN(OP,SCALAR_IMPL,EXPR_IMPL)                \
    template <typename OtherT>                                          \
    Array const &                                                       \
    operator OP(Expression<OtherT> const & expr) const {                \
        NDARRAY_ASSERT(expr.getShape() == this->getShape());            \
        EXPR_IMPL(OP);                                                  \
        return *this;                                                   \
    }                                                                   \
    template <typename ScalarT>                                         \
    typename boost::enable_if<boost::is_convertible<ScalarT,Element>, Array const &>::type \
    operator OP(ScalarT const & scalar) const {                         \
        SCALAR_IMPL(OP);                                                \
        return *this;                                                   \
    }
#else // DOXYGEN
#define NDARRAY_GENERAL_ASSIGN(OP,SCALAR_IMPL,EXPR_IMPL)                \
    template <typename OtherT>                                          \
    Array const &                                                       \
    operator OP(Expression<OtherT> const & expr) const;                 \
    template <typename ScalarT>                                         \
    Array const &                                                       \
    operator OP(ScalarT const & scalar) const;
#endif

#define NDARRAY_BASIC_ASSIGN_SCALAR(OP)         \
    std::fill(this->begin(),this->end(),scalar)

#define NDARRAY_BASIC_ASSIGN_EXPR(OP)                   \
    std::copy(expr.begin(),expr.end(),this->begin())

#define NDARRAY_AUGMENTED_ASSIGN_SCALAR(OP)                             \
    Iterator const i_end = this->end();                                 \
    for (Iterator i = this->begin(); i != i_end; ++i) (*i) OP scalar

#define NDARRAY_AUGMENTED_ASSIGN_EXPR(OP)                               \
    Iterator const i_end = this->end();                                 \
    typename OtherT::Iterator j = expr.begin();                         \
    for (Iterator i = this->begin(); i != i_end; ++i, ++j) (*i) OP (*j)

#define NDARRAY_BASIC_ASSIGN                                            \
    NDARRAY_GENERAL_ASSIGN(=,NDARRAY_BASIC_ASSIGN_SCALAR,NDARRAY_BASIC_ASSIGN_EXPR)

#define NDARRAY_AUGMENTED_ASSIGN(OP)                                    \
    NDARRAY_GENERAL_ASSIGN(OP,NDARRAY_AUGMENTED_ASSIGN_SCALAR,NDARRAY_AUGMENTED_ASSIGN_EXPR)

namespace ndarray {

/**
 *  @class Array
 *  \brief Multidimensional array
 *
 *  \ingroup MainGroup
 *
 *  Array is the workhorse multidimensional array class for ndarray.  It generally mimics a
 *  container of nested lower-dimensional arrays (though it does not completely conform to the STL
 *  Container concept).  Array is designed to allow extremely lightweight nested iterators
 *  an nested array extraction, and provides an optional template parameter that enables
 *  additional optimizations for contiguous row-major arrays.
 *
 *  Array objects provide a shallow, reference-counted copy constructor and deep-copying assignment
 *  operators (along with a special syntax for shallow assignment; see shallow()).
 *  Construction of new arrays is handled by the allocate() and external() functions.  
 *
 *  Because memory is shared between arrays, Array supports constness-propogation semantics more similar
 *  to smart pointer objects such as auto_ptr or shared_ptr than standard library containers.  The Array
 *  template can be instantiated with a const data type, which makes it essentially an immutable container;
 *  but the constness of an Array object or reference does not affect the mutability of the underlying data.
 *  As a result:
 *  \li <b><tt>Array<int,3></tt></b> allows both deep and shallow assignment.
 *  \li <b><tt>Array<int,3> const</tt></b> allows deep assignment but not shallow assignment.
 *  \li <b><tt>Array<int const,3></tt></b> allows shallow assignment but not deep assignment.
 *  \li <b><tt>Array<int const,3> const</tt></b> allows neither deep or shallow assignment.
 *
 *  This has a few somewhat suprising implications:
 *  \li Assignment and augmented assignment operators for Array are const member functions.
 *  \li All element-access member functions are const.
 *  \li Array output parameters will typically be passed by const reference unless shallow assignment is
 *      desired.
 *
 *  \tparam T Data type of array elements.
 *  \tparam N Number of dimensions.
 *  \tparam C Number of guaranteed row-major contiguous dimensions, starting with the last.  Defaults to 0.
 *
 */
template <typename T, int N, int C>
class Array : public detail::ArrayImpl< Array<T,N,C> > {
    typedef detail::ArrayImpl<Array> Super;
    typedef typename detail::Core<typename boost::remove_const<T>::type,N> Core;
    typedef typename Core::ConstPtr CorePtr;
public:
    /// \brief Data type of array elements.
    typedef typename detail::ArrayTraits<Array>::Element Element;
    /// \brief Number of dimensions (boost::mpl::int_).
    typedef typename detail::ArrayTraits<Array>::ND ND;
    /// \brief Number of guaranteed row-major contiguous dimensions, counted from the end (boost::mpl::int_).
    typedef typename detail::ArrayTraits<Array>::RMC RMC;
    /// \brief Nested array or element iterator.
    typedef typename detail::ArrayTraits<Array>::Iterator Iterator; 
    /// \brief Nested array or element reference.
    typedef typename detail::ArrayTraits<Array>::Reference Reference;
    /// \brief Nested array or element value type.
    typedef typename detail::ArrayTraits<Array>::Value Value;
    /// \brief Shared pointer that manages memory.
    typedef boost::shared_ptr<Element> Owner;
    /// \brief Vector type for N-dimensional indices.
    typedef Vector<int,N> Index;

    /** 
     *  \brief Default constructor. 
     *
     *  Creates an empty array with zero dimensions and null memory.
     */
    Array() : _data(0), _core(Core::create()) {}

    /** 
     *  \brief Converting copy constructor. 
     *
     *  Implicit conversion is allowed for non-const -> const and for 
     *  more guaranteed RMC -> less guaranteed RMC (see \ref overview).
     */
    template <typename T1, int C1>
    Array(
        Array<T1,N,C1> const & other
#ifndef DOXYGEN
        , typename boost::enable_if_c<((C1>=C) && boost::is_convertible<T1*,T*>::value),void*>::type=0
#endif
    ) : _data(other._data), _core(other._core) {}
    
    /// \brief Return a raw pointer to the beginning of the array.
    Element * getData() const { return this->_data; }

    /// \brief Return the shared_ptr that manages the lifetime of the array data.
    Owner getOwner() const { return this->_core->getOwner(); }

    /// \brief Return the size of a specific dimension.
    template <int P> int getSize() const {
        return detail::getDimension<P>(*this->_core).getSize();
    }

    /// \brief Return the stride in a specific dimension.
    template <int P> int getStride() const {
        return detail::getDimension<P>(*this->_core).getStride();
    }

    /// \brief Return a Vector of the sizes of all dimensions.
    Index getShape() const { Index r; this->_core->fillShape(r); return r; }

    /// \brief Return a Vector of the strides of all dimensions.
    Index getStrides() const { Index r; this->_core->fillStrides(r); return r; }

    /// \briref Return the total number of elements in the array.
    int getNumElements() const {
        return this->_core->getNumElements();
    }

    /// \internal \brief A template metafunction class to determine the result of a view indexing operation.
    template <
        typename ViewT, 
        typename SeqT = typename detail::ViewNormalizer<N,typename ViewT::Sequence>::Output
        >
    struct ResultOf {
        typedef T Element;
        typedef typename detail::ViewTraits<N,C,SeqT>::ND ND;
        typedef typename detail::ViewTraits<N,C,SeqT>::RMC RMC;
        typedef Array<Element,ND::value,RMC::value> Type;
    };

    /// \brief Return a general view into this array (see \ref tutorial).
    template <typename ViewT>
    typename ResultOf<ViewT>::Type
    operator[](ViewT const & def) const {
        return detail::buildView(*this,def._seq);
    }

    /// \brief Return a single subarray (for ND > 1) or element for (ND==1).
    using detail::ArrayImpl<Array>::operator[];

    /// \brief Deep assignment operator.
    Array const & operator=(Array const & other) const {
        if (&other != this) {
            NDARRAY_ASSERT(other.getShape() == this->getShape());
            std::copy(other.begin(),other.end(),this->begin());
        }
        return *this;
    }

    NDARRAY_BASIC_ASSIGN
    NDARRAY_AUGMENTED_ASSIGN(+=)
    NDARRAY_AUGMENTED_ASSIGN(-=)
    NDARRAY_AUGMENTED_ASSIGN(*=)
    NDARRAY_AUGMENTED_ASSIGN(/=)
    NDARRAY_AUGMENTED_ASSIGN(%=)
    NDARRAY_AUGMENTED_ASSIGN(^=)
    NDARRAY_AUGMENTED_ASSIGN(&=)
    NDARRAY_AUGMENTED_ASSIGN(|=)
    NDARRAY_AUGMENTED_ASSIGN(<<=)
    NDARRAY_AUGMENTED_ASSIGN(>>=)

private:
    template <typename T1, int N1, int C1> friend class Array;
    template <typename ArrayT1> friend class detail::ArrayAccess;
    template <typename ArrayT1, typename SeqT1> friend class detail::ViewBuilder;
    template <int N1, typename Allocator1> friend class detail::AllocationInitializer;
    template <typename T1, int N1> friend class detail::ExternalInitializer;

    /// \internal \brief Construct an Array from a pointer and Core.
    Array(Element * data, CorePtr const & core) : _data(data), _core(core) {}

    Element * _data;
    CorePtr _core;
};

} // namespace ndarray

#undef NDARRAY_BASIC_ASSIGN
#undef NDARRAY_BASIC_ASSIGN_SCALAR
#undef NDARRAY_BASIC_ASSIGN_EXPR
#undef NDARRAY_AUGMENTED_ASSIGN
#undef NDARRAY_AUGMENTED_ASSIGN_SCALAR
#undef NDARRAY_AUGMENTED_ASSIGN_EXPR
#undef NDARRAY_GENERAL_ASSIGN

#endif // !NDARRAY_Array_hpp_INCLUDED
