#ifndef NDARRAY_casts_hpp_INCLUDED
#define NDARRAY_casts_hpp_INCLUDED

/** 
 *  @file ndarray/casts.hpp
 *
 *  \brief Specialized casts for Array.
 */

#include "ndarray/Array.hpp"
#include <boost/type_traits/add_const.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/mpl/comparison.hpp>
#include <boost/static_assert.hpp>

namespace ndarray {

namespace detail {

template <typename Original, typename Casted>
class ReinterpretDeleter {
    boost::shared_ptr<Original> _original;
public:
    void operator()(Casted * p) { _original.reset(); }
    explicit ReinterpretDeleter(boost::shared_ptr<Original> const & original) : _original(original) {}

    static boost::shared_ptr<Casted> makeOwner(boost::shared_ptr<Original> const & original) {
        return boost::shared_ptr<Casted>(
            reinterpret_cast<Casted*>(original.get()), 
            ReinterpretDeleter(original)
        );
    }
};

template <typename ArrayT>
struct ComplexExtractor {
    typedef typename detail::ArrayTraits<ArrayT>::Element Element;
    typedef typename boost::remove_const<Element>::type NonConst;
    typedef typename boost::is_const<Element>::type IsConst;
    BOOST_STATIC_ASSERT( boost::is_complex<NonConst>::value );
    typedef typename NonConst::value_type NonConstValue;
    typedef typename boost::add_const<NonConstValue>::type ConstValue;
    typedef typename boost::mpl::if_<IsConst, ConstValue, NonConstValue>::type Value;
    typedef typename detail::ArrayTraits<ArrayT>::ND ND;
    typedef Array<Value,ND::value,0> View;
    typedef Vector<int,ND::value> Index;
    typedef ReinterpretDeleter<Element,Value> Deleter;

    static inline View apply(ArrayT const & array, int offset) {
        return View(
            external(
                reinterpret_cast<Value*>(array.getData()) + offset,
                array.getShape(),
                array.getStrides() * 2,
                Deleter::makeOwner(array.getOwner())
            )
        );
    }
};

} // namespace detail

/// \addtogroup MainGroup
/// @{

/**
 *  Convert an Array with a const data type to an array
 *  with a non-const data type.
 */
template <typename T, typename ArrayT>
typename detail::ArrayTraits<ArrayT>::template RebindT<T>::Other
const_array_cast(ArrayT const & array) {
    typedef typename detail::ArrayTraits<ArrayT>::template RebindT<T>::Other Output;
    return detail::ArrayAccess<Output>::construct(
        const_cast<typename Output::Element*>(array.getData()),
        array
    );
}

/**
 *  Convert an Array to a type with more guaranteed
 *  row-major-contiguous dimensions with no checking.
 */
template <int C, typename ArrayT>
typename detail::ArrayTraits<ArrayT>::template RebindC<C>::Other
static_dimension_cast(ArrayT const & array) {
    typedef typename detail::ArrayTraits<ArrayT>::template RebindC<C>::Other Output;
    return detail::ArrayAccess<Output>::construct(
        const_cast<typename Output::Element*>(array.getData()),
        array
    );
}

/**
 *  Convert an Array to a type with more guaranteed
 *  row-major-contiguous dimensions, if the strides
 *  of the array match the desired number of RMC
 *  dimensions.  If the cast fails, an empty Array
 *  is returned.
 */
template <int C, typename ArrayT>
typename detail::ArrayTraits<ArrayT>::template RebindC<C>::Other
dynamic_dimension_cast(ArrayT const & array) {
    typedef typename detail::ArrayTraits<ArrayT>::template RebindC<C>::Other Output;
    const static int N = Output::ND::value;
    Vector<int,N> shape = array.getShape();
    Vector<int,N> strides = array.getStrides();
    int n = 1;
    for (int i=1; i<=C; ++i) {
        if (strides[N-i] != n) return Output();
        n *= shape[N-i];
    }
    return static_dimension_cast<C>(array);
}

template <typename ArrayT>
typename detail::ComplexExtractor<ArrayT>::View
getReal(ArrayT const & array) {
    return detail::ComplexExtractor<ArrayT>::apply(array, 0);
}

template <typename ArrayT>
typename detail::ComplexExtractor<ArrayT>::View
getImag(ArrayT const & array) {
    return detail::ComplexExtractor<ArrayT>::apply(array, 1);
}

/// @}

}

#endif // !NDARRAY_casts_hpp_INCLUDED
