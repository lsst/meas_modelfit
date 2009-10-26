#ifndef NDARRAY_DETAIL_ArrayImpl_hpp_INCLUDED
#define NDARRAY_DETAIL_ArrayImpl_hpp_INCLUDED

/** 
 *  @file ndarray/detail/ArrayImpl.hpp
 *
 *  \brief Definitions for ArrayImpl.
 */

#include "ndarray/Expression.hpp"
#include "ndarray/detail/NestedIterator.hpp"
#include "ndarray/detail/StridedIterator.hpp"
#include "ndarray/detail/ArrayAccess.hpp"

namespace ndarray {
namespace detail {

/**
 *  \internal @class ArrayImpl
 *  \brief CRTP implementation for Array.
 *
 *  \ingroup InternalGroup
 *
 *  Implements member functions that need specialization
 *  for 1D arrays.
 */
template <typename Derived, int N, int C>
class ArrayImpl : public Expression<Derived> {
public:
    typedef typename ArrayTraits<Derived>::Element Element;
    typedef typename ArrayTraits<Derived>::Reference Reference;
    typedef typename ArrayTraits<Derived>::Iterator Iterator;

    /// \brief Return a single subarray.
    Reference operator[](int n) const {
        return ArrayAccess<Reference>::construct(
            this->getSelf().getData() + n * this->getSelf().template getStride<0>(), 
            this->getSelf()
        );
    }

    /// \brief Return an Iterator to the beginning of the array.
    Iterator begin() const {
        return Iterator(
            ArrayAccess<Reference>::construct(this->getSelf().getData(), this->getSelf()), 
            this->getSelf().template getStride<0>()
        );
    }

    /// \brief Return an Iterator to one past the end of the array.
    Iterator end() const {
        return Iterator(
            ArrayAccess<Reference>::construct(
                this->getSelf().getData() 
                + this->getSelf().template getSize<0>() * this->getSelf().template getStride<0>(), 
                this->getSelf()
            ),
            this->getSelf().template getStride<0>()
        );
    }

protected:
    ArrayImpl() {}
};

/**
 *  \internal
 *  \brief CRTP implementation for Array.
 *
 *  \ingroup InternalGroup
 *
 *  Implements member functions that need specialization
 *  for 1D arrays (1D contiguous specialization).
 */
template <typename Derived>
class ArrayImpl<Derived,1,1> : public Expression<Derived> {
public:
    typedef typename ArrayTraits<Derived>::Element Element;
    typedef typename ArrayTraits<Derived>::Reference Reference;
    typedef typename ArrayTraits<Derived>::Iterator Iterator;

    /// \brief Return a single element.
    Reference operator[](int n) const {
        return *(this->getSelf().getData() + n);
    }
    
    /// \brief Return an Iterator to the beginning of the array.
    Iterator begin() const {
        return this->getSelf().getData(); 
    }

    /// \brief Return an Iterator to one past the end of the array.
    Iterator end() const {
        return this->getSelf().getData() + this->getSelf().template getSize<0>();
    }

protected:
    ArrayImpl() {}
};


/**
 *  \internal
 *  \brief CRTP implementation for Array.
 *
 *  \ingroup InternalGroup
 *
 *  Implements member functions that need specialization
 *  for 1D arrays (1D noncontiguous specialization).
 */
template <typename Derived>
class ArrayImpl<Derived,1,0> : public Expression<Derived> {
public:
    typedef typename ArrayTraits<Derived>::Element Element;
    typedef typename ArrayTraits<Derived>::Reference Reference;
    typedef typename ArrayTraits<Derived>::Iterator Iterator;

    /// \brief Return a single element.
    Reference operator[](int n) const {
        return *(this->getSelf().getData() + n * this->getSelf().template getStride<0>());
    }

    /// \brief Return an Iterator to the beginning of the array.
    Iterator begin() const {
        return Iterator(
            this->getSelf().getData(),
            this->getSelf().template getStride<0>()
        ); 
    }

    /// \brief Return an Iterator to one past the end of the array.
    Iterator end() const {
        return Iterator(
            this->getSelf().getData() 
            + this->getSelf().template getSize<0>() * this->getSelf().template getStride<0>(),
            this->getSelf().template getStride<0>()
        );
    }

protected:
    ArrayImpl() {}
};

} // namespace ndarray::detail
} // namespace ndarray

#endif // !NDARRAY_DETAIL_ArrayImpl_hpp_INCLUDED
