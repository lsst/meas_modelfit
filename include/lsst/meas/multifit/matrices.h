// -*- lsst-c++ -*-
/**
 * @file
 *
 * Collection of utility functions for Eigen and ndarray
 */
#ifndef LSST_MEAS_MULTIFIT_MATRICES_H
#define LSST_MEAS_MULTIFIT_MATRICES_H

#include <Eigen/Core>
#include <ndarray.hpp>

#include "lsst/meas/multifit/core.h"

namespace lsst {
namespace meas {
namespace multifit {

#if 0 
// TODO: Reimplement this when Eigen::Map has a constructor that takes strides.
/**
 * Get an 2d Eigen::Map over a 2d ndarray
 */
inline MatrixMapBlock getMatrixView(
    ndarray::Array<Pixel const,2,1> const & array
) {
    MatrixMap map(array.getData(), array.getStride<0>(), array.getSize<0>());
    return MatrixMapBlock(map, 0, 0, array.getSize<1>(), array.getSize<0>());
}
#endif

/**
 * Get an 1d EigenMap over a 1d ndarray
 */
inline VectorMap getVectorView(ndarray::Array<Pixel const,1,1> const & array) {
    return VectorMap(array.getData(), array.getSize<0>());
}

/**
 * Get a 2d Eigen::Map over a 3d ndarray
 *
 * The innermost dimensions are flattened, such that the resulting map has
 * (array.getSize<2>() * array.getSize<1>()) rows, and array.getSize<0> columns
 */
template <typename T>
inline Eigen::Map< Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> >
getCompressedMatrixView(ndarray::Array<T,3,3> const & array) {
    return Eigen::Map< Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> >(
        array.getData(),
        array.template getSize<2>() * array.template getSize<1>(),
        array.template getSize<0>()
    );
}

/**
 * Get a 1d Eigen::Map over a 2d ndarray
 *
 * The dimensions of the ndarray are flattened, such that the resulting map has
 * (array.getSize<1>()*array.getSize<0>() rows, and 1 column.
 */
template <typename T>
inline Eigen::Map< Eigen::Matrix<T,Eigen::Dynamic,1> >
getCompressedVectorView(ndarray::Array<T,2,2> const & array) {
    return Eigen::Map< Eigen::Matrix<T,Eigen::Dynamic,1> >(
        array.getData(),
        array.template getSize<1>() * array.template getSize<0>()
    );
}

/**
 * Utility wrapper for obtaining a view of a 2d ndarray as specified by a BoxI
 *
 * @return a view of input which spans the rows [bbox.getMinY(), 
 *      bbox.getMaxY()] and columns [bbox.getMinx(), bbox.getMaxX()]
 */
template <typename T, int C>
inline ndarray::Array<T, 2, ((C>=1) ? 1:0)> window(
    ndarray::Array<T, 2, C> const & input, 
    lsst::afw::geom::BoxI const & box
) {
    return input[ndarray::view
        (box.getMinY(), box.getMaxY()+1)
        (box.getMinX(), box.getMaxX()+1)
    ];
}

/**
 * Utility wrapper for obtaining a view of a 3d ndarray as specified by a BoxI
 *
 * @return a view of input which spans the includes the entire range of the
 *      outer most dimension, rows [bbox.getMinY(), 
 *      bbox.getMaxY()] and columns [bbox.getMinx(), bbox.getMaxX()]
 */
template <typename T, int C>
inline ndarray::Array<T, 3, ((C>=1) ? 1:0)> window(
    ndarray::Array<T, 3, C> const & input, 
    lsst::afw::geom::BoxI const & box
) {
    return input[ndarray::view
        ()
        (box.getMinY(), box.getMaxY()+1)
        (box.getMinX(), box.getMaxX()+1)
    ];
}

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_MATRICES_H
