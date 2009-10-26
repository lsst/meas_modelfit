#ifndef NDARRAY_DETAIL_NestedIterator_hpp_INCLUDED
#define NDARRAY_DETAIL_NestedIterator_hpp_INCLUDED

/** 
 *  @file ndarray/detail/NestedIterator.hpp
 *
 *  \brief Definition of NestedIterator.
 */

#include <boost/iterator/iterator_facade.hpp>
#include "ndarray/detail/ArrayTraits.hpp"
#include "ndarray/shallow.hpp"

namespace ndarray {
namespace detail {

/**
 *  \internal \brief Nested-array iterator class for Array with ND > 1.
 *
 *  While this iterator makes use of boost::iterator_facade, it
 *  reimplements the actual dereferencing operations (operator*,
 *  operator->) to return <b><tt>Reference const &</tt></b> and
 *  <b><tt>Reference const *</tt></b>, even though these are
 *  only convertible to the <b><tt>reference</tt></b> and
 *  <b><tt>pointer</tt></b> types associated with the iterator,
 *  not the types themselves.
 *
 *  \ingroup InternalGroup
 */
template <typename ArrayT>
class NestedIterator : public boost::iterator_facade<
    NestedIterator<ArrayT>,
    typename ArrayTraits<ArrayT>::Value,
    boost::random_access_traversal_tag,
    typename ArrayTraits<ArrayT>::Reference
    >
{
public:
    typedef typename ArrayTraits<ArrayT>::Value Value;
    typedef typename ArrayTraits<ArrayT>::Reference Reference;
    
    Reference operator[](int n) const {
        Reference r(_ref);
        ArrayAccess<Reference>::getArrayData(r) += n * _stride;
        return r;
    }

    Reference const & operator*() const { return _ref; }

    Reference const * operator->() { return &_ref; }

    NestedIterator() : _ref(), _stride(0) {}

    NestedIterator(Reference const & ref, int stride) : _ref(ref), _stride(stride) {}

    NestedIterator(NestedIterator const & other) : _ref(other._ref), _stride(other._stride) {}

    template <typename OtherT>
    NestedIterator(NestedIterator<OtherT> const & other) : _ref(other._ref), _stride(other._stride) {}

    NestedIterator & operator=(NestedIterator const & other) {
        if (&other != this) {
            ndarray::shallow(_ref) = other._ref;
            _stride = other._stride;
        }
        return *this;
    }

    template <typename OtherT>
    NestedIterator & operator=(NestedIterator<OtherT> const & other) {
        ndarray::shallow(_ref) = other._ref;
        _stride = other._stride;
        return *this;
    }

private:

    friend class boost::iterator_core_access;

    template <typename OtherT> friend class NestedIterator;

    Reference const & dereference() const { return _ref; }

    void increment() { ArrayAccess<Reference>::getArrayData(_ref) += _stride; }
    void decrement() { ArrayAccess<Reference>::getArrayData(_ref) -= _stride; }
    void advance(int n) { ArrayAccess<Reference>::getArrayData(_ref) += _stride * n; }

    template <typename OtherT>
    int distance_to(NestedIterator<OtherT> const & other) const {
        return std::distance(_ref.getData(), other._ref.getData()) / _stride; 
    }

    template <typename OtherT>
    bool equal(NestedIterator<OtherT> const & other) const {
        return _ref.getData() == other._ref.getData();
    }

    Reference _ref;
    int _stride;
};

} // namespace ndarray::detail
} // namespace ndarray

#endif // !NDARRAY_DETAIL_NestedIterator_hpp_INCLUDED
