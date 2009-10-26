#ifndef NDARRAY_shallow_hpp_INCLUDED
#define NDARRAY_shallow_hpp_INCLUDED

/** 
 *  \file ndarray/shallow.hpp \brief Shallow assignment and comparison for Array.
 */

#include "ndarray/Array.hpp"

namespace ndarray {
namespace detail {

/** 
 *  \internal
 *  \ingroup InternalGroup
 *  \brief An expression template that provides shallow comparison for Array.
 */
template <typename ArrayT>
class Shallow<ArrayT const> {
public:

    explicit Shallow(ArrayT const & array) : _array(const_cast<ArrayT&>(array)) {}

    template <typename OtherT>
    bool operator==(Shallow<OtherT> const & other) const {
        return _array.getData() == other._array.getData()
            && _array.getShape() == other._array.getShape()
            && _array.getStrides() == other._array.getStrides();
    }

    template <typename OtherT>
    bool operator!=(Shallow<OtherT> const & other) const {
        return !this->operator==(other);
    }

protected:
    template <typename OtherT> friend class Shallow;
    ArrayT & _array;
private:
    void operator=(Shallow const & other) {}
};

/** 
 *  \internal
 *  \ingroup InternalGroup
 *  \brief An expression template that provides shallow comparison and assignment for Array.
 */
template <typename ArrayT> 
class Shallow : public Shallow<ArrayT const> {
    typedef Shallow<ArrayT const> Super;
public:

    explicit Shallow(ArrayT & array) : Super(array) {}

    ArrayT & operator=(ArrayT const & other) {
        if (&other != &this->_array) ArrayAccess<ArrayT>::assign(this->_array,other);
        return this->_array;
    }

    void reset() {
        this->operator=(ArrayT());
    }

};

} // namespace ndarray::detail

/** 
 *  \brief Create a temporary expression that allows for shallow comparison and assignment of arrays.
 *
 *  The shallow function provides a syntax for shallow comparisons and assignments:
 *  \code
 *  Array<int,2> a = allocate(makeVector(4,5));
 *  Array<int,2> b = allocate(makeVector(4,5));
 *  assert(shallow(a) != shallow(b));
 *  shallow(b) = a;
 *  assert(shallow(a) == shallow(b));
 *  \endcode
 *  Two arrays are considered equal in a shallow sense if they have identical shape and strides
 *  and refer to the same block of memory.
 *
 *  \ingroup MainGroup
 */
template <typename ArrayT>
inline detail::Shallow<ArrayT>
shallow(ArrayT & array) {
    return detail::Shallow<ArrayT>(array);
}

/** 
 *  \brief Create a temporary expression that allows for shallow comparison of arrays.
 *
 *  This overload, for const Array references, supports comparison but not assignment.
 *
 *  \ingroup MainGroup
 *  \sa shallow
 */
template <typename ArrayT>
inline detail::Shallow<ArrayT const> 
shallow(ArrayT const & array) {
    return detail::Shallow<ArrayT const>(array);
}

} // namespace ndarray

#endif // !NDARRAY_shallow_hpp_INCLUDED
