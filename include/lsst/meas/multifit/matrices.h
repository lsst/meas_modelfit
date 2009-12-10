#ifndef LSST_MEAS_MULTIFIT_MATRICES_H
#define LSST_MEAS_MULTIFIT_MATRICES_H

#include <Eigen/Core>
#include <ndarray.hpp>

#include "lsst/meas/multifit/core.h"

namespace lsst {
namespace meas {
namespace multifit {

typedef Eigen::aligned_allocator<char> Allocator;

typename Eigen::Map< Eigen::Matrix<Pixel,Eigen::Dynamic,Eigen::Dynamic> > MatrixMap;
typename Eigen::Map< Eigen::Matrix<Pixel,1,Eigen::Dynamic> > VectorMap;

inline MatrixMap getCompressedMatrixView(ndarray::Array<Pixel,2,2> const & array) {
    return MatrixMap(array.getData(), array.getSize<1>(), array.getSize<0>());
}

inline VectorMap getCompressedVectorView(ndarray::Array<Pixel,1,1> const & array) {
    return VectorMap(array.getData(), array.getSize<0>());
}

template <typename T>
inline Eigen::Map< Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> >
getMatrixView(ndarray::Array<T,3,3> const & array) {
    return Eigen::Map< Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> >(
        array.getData(),
        array.template getSize<2>() * array.template getSize<1>(),
        array.template getSize<0>()
    );
}

template <typename T>
inline Eigen::Map< Eigen::Matrix<T,Eigen::Dynamic,1> >
getVectorView(ndarray::Array<T,2,2> const & array) {
    return Eigen::Map< Eigen::Matrix<T,Eigen::Dynamic,1> >(
        array.getData(),
        array.template getSize<1>() * array.template getSize<0>()
    );
}

template <typename T, int C>
inline ndarray::Array<T,2,((C>=1) ? 1:0)> window(
    ndarray::Array<T,2,C> const & input, 
    lsst::afw::image::BBox const & box
) {
    return input[ndarray::view
        (box.getY0(),box.getY1()+1)
        (box.getX0(),box.getX1()+1)
    ];
}

template <typename T, int C>
inline ndarray::Array<T,3,((C>=1) ? 1:0)> window(
    ndarray::Array<T,3,C> const & input, 
    lsst::afw::image::BBox const & box
) {
    return input[ndarray::view
        ()
        (box.getY0(),box.getY1()+1)
        (box.getX0(),box.getX1()+1)
    ];
}

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_MATRICES_H
