#ifndef LSST_MEAS_MULTIFIT_MATRICES_H
#define LSST_MEAS_MULTIFIT_MATRICES_H

#include <Eigen/Core>
#include <ndarray.hpp>

#include "lsst/afw/geom/Box.h"
#include "lsst/meas/multifit/core.h"

namespace lsst {
namespace meas {
namespace multifit {

typedef Eigen::aligned_allocator<char> Allocator;

typedef Eigen::Map<Eigen::Matrix<Pixel, Eigen::Dynamic, Eigen::Dynamic> > MatrixMap;
typedef Eigen::Block<MatrixMap> MatrixMapBlock;
typedef Eigen::Map<Eigen::Matrix<Pixel, 1, Eigen::Dynamic> > VectorMap;


inline MatrixMapBlock getMatrixView(ndarray::Array<Pixel const,2,1> const & array) {
    MatrixMap map(array.getData(), array.getStride<0>(), array.getSize<0>());
    return MatrixMapBlock(map, 0, 0, array.getSize<1>(), array.getSize<0>());
}

inline VectorMap getVectorView(ndarray::Array<Pixel const,1,1> const & array) {
    return VectorMap(array.getData(), array.getSize<0>());
}

template <typename T>
inline Eigen::Map< Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> >
getCompressedMatrixView(ndarray::Array<T,3,3> const & array) {
    return Eigen::Map< Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> >(
        array.getData(),
        array.template getSize<2>() * array.template getSize<1>(),
        array.template getSize<0>()
    );
}

template <typename T>
inline Eigen::Map< Eigen::Matrix<T,Eigen::Dynamic,1> >
getCompressedVectorView(ndarray::Array<T,2,2> const & array) {
    return Eigen::Map< Eigen::Matrix<T,Eigen::Dynamic,1> >(
        array.getData(),
        array.template getSize<1>() * array.template getSize<0>()
    );
}

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
