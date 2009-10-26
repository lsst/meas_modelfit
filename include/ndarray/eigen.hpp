#ifndef NDARRAY_eigen_hpp_INCLUDED
#define NDARRAY_eigen_hpp_INCLUDED

/**
 *  @file ndarray/eigen.hpp
 *  \brief Conversions between ndarray::Array and Eigen objects.
 *
 *  \note This file is not included by the main "ndarray.hpp" header file.
 */

#include "ndarray.hpp"
#include <Eigen/Core>

namespace ndarray {
namespace detail {

/**
 *  \internal \ingroup InternalGroup
 *  \brief Helper class for Array-to-Eigen conversions.
 */
template <typename Scalar, int Rows, int Cols>
struct ArrayAsEigenMap {
    typedef Eigen::Matrix<Scalar,Rows,Cols,(Eigen::RowMajor | Eigen::DontAlign)> Matrix;
    typedef Eigen::Map<Matrix> Map;
    
    static Map apply(Array<Scalar,2,2> const & array) {
        return Map(array.getData(), array.template getSize<0>(), array.template getSize<1>());
    }

    static Map const apply(Array<Scalar const,2,2> const & array) {
        return Map(array.getData(), array.template getSize<0>(), array.template getSize<1>());
    }

    static Map apply(Array<Scalar,1,1> const & array) {
        NDARRAY_ASSERT(false);
        return Map(0,0,0);
    }

    static Map const apply(Array<Scalar const,1,1> const & array) {
        NDARRAY_ASSERT(false);
        return Map(0,0,0);
    }
    
};

/// \cond SPECIALIZATIONS

/**
 *  \internal \ingroup InternalGroup
 *  \brief Helper class for Array-to-Eigen conversions (Eigen::VectorN specialization).
 */
template <typename Scalar, int Rows>
struct ArrayAsEigenMap<Scalar,Rows,1> {
    typedef Eigen::Matrix<Scalar,Rows,1,(Eigen::RowMajor | Eigen::DontAlign)> Matrix;
    typedef Eigen::Map<Matrix> Map;
    
    static Map apply(Array<Scalar,2,2> const & array) {
        return Map(array.getData(), array.template getSize<0>(), array.template getSize<1>());
    }

    static Map const apply(Array<Scalar const,2,2> const & array) {
        return Map(array.getData(), array.template getSize<0>(), array.template getSize<1>());
    }

    static Map apply(Array<Scalar,1,1> const & array) {
        return Map(array.getData(), array.size());
    }

    static Map const apply(Array<Scalar const,1,1> const & array) {
        return Map(array.getData(), array.size());
    }
    
};

/**
 *  \internal \ingroup InternalGroup
 *  \brief Helper class for Array-to-Eigen conversions (Eigen::RowVectorN specialization).
 */
template <typename Scalar, int Cols>
struct ArrayAsEigenMap<Scalar,1,Cols> {
    typedef Eigen::Matrix<Scalar,1,Cols,(Eigen::RowMajor | Eigen::DontAlign)> Matrix;
    typedef Eigen::Map<Matrix> Map;
    
    static Map apply(Array<Scalar,2,2> const & array) {
        return Map(array.getData(), array.template getSize<0>(), array.template getSize<1>());
    }

    static Map const apply(Array<Scalar const,2,2> const & array) {
        return Map(array.getData(), array.template getSize<0>(), array.template getSize<1>());
    }

    static Map apply(Array<Scalar,1,1> const & array) {
        return Map(array.getData(), array.size());
    }

    static Map const apply(Array<Scalar const,1,1> const & array) {
        return Map(array.getData(), array.size());
    }
    
};

/// \endcond

} // namespace ndarray::eigen::detail

#ifndef DOXYGEN
template <typename Matrix>
typename boost::enable_if_c< 
    (Matrix::Flags & Eigen::DirectAccessBit) && (Matrix::Flags & Eigen::LinearAccessBit), 
    Array<typename boost::mpl::if_<boost::is_const<Matrix>,
                                   typename Matrix::Scalar const, 
                                   typename Matrix::Scalar
                                   >::type,
          1>
>::type
#else
/**
 *  \ingroup MainGroup
 *  \brief Create a 1D Array view into an Eigen object.
 *
 *  The created Array does not own a reference to its data, so the user is responsible for 
 *  ensuring the memory remains valid for the lifetime of the array.
 */
Array<typename Matrix::Scalar,1>
#endif
viewVectorAsArray(Matrix & matrix) {
    return external(
        const_cast<typename Matrix::Scalar*>(matrix.data()),
        makeVector(matrix.size()),
        makeVector(1)
    );
}

#ifndef DOXYGEN
template <typename Matrix>
typename boost::enable_if_c< 
    Matrix::Flags & Eigen::DirectAccessBit, 
    Array<typename boost::mpl::if_<boost::is_const<Matrix>,
                                   typename Matrix::Scalar const, 
                                   typename Matrix::Scalar
                                   >::type,
          2>
>::type
#else
/**
 *  \ingroup MainGroup
 *  \brief Create a 2D Array view into an Eigen object.
 *
 *  The created Array does not own a reference to its data, so the user is responsible for 
 *  ensuring the memory remains valid for the lifetime of the array.
 */
Array<typename Matrix::Scalar,2>
#endif
viewMatrixAsArray(Matrix & matrix) {
    if (Matrix::Flags & Eigen::RowMajorBit) {
        return external(
            const_cast<typename Matrix::Scalar*>(matrix.data()),
            makeVector(matrix.rows(),matrix.cols()),
            makeVector(matrix.stride(),1)
        );
    } else {
        return external(
            const_cast<typename Matrix::Scalar*>(matrix.data()),
            makeVector(matrix.rows(),matrix.cols()),
            makeVector(1,matrix.stride())
        );        
    }
}

/**
 *  \ingroup MainGroup
 *  \brief Create an Eigen::Map to an Array.
 *
 *  The first template parameter specifies an Eigen type
 *  compatible with the returned Map object, and must
 *  be explicitly provided.
 */
template <typename Matrix, int N>
typename detail::ArrayAsEigenMap<typename Matrix::Scalar,
                                 Matrix::RowsAtCompileTime,
                                 Matrix::ColsAtCompileTime>::Map
viewArrayAs(Array<typename Matrix::Scalar,N,N> const & array) {
    return detail::ArrayAsEigenMap<typename Matrix::Scalar,
        Matrix::RowsAtCompileTime,
        Matrix::ColsAtCompileTime>::apply(array);
}

/// \cond SPECIALIZATIONS
/**
 *  \ingroup MainGroup
 *  \brief Create an Eigen::Map to an Array (const specialization).
 *
 *  The first template parameter specifies an Eigen type
 *  compatible with the returned Map object, and must
 *  be explicitly provided.
 */
template <typename Matrix, int N>
typename detail::ArrayAsEigenMap<typename Matrix::Scalar,
                                 Matrix::RowsAtCompileTime,
                                 Matrix::ColsAtCompileTime>::Map const
viewArrayAs(Array<typename Matrix::Scalar const,N,N> const & array) {
    return detail::ArrayAsEigenMap<typename Matrix::Scalar,
        Matrix::RowsAtCompileTime,
        Matrix::ColsAtCompileTime>::apply(array);
}
/// \endcond

} // namespace ndarray

#endif // !NDARRAY_eigen_hpp_INCLUDED
