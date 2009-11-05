#ifndef LSST_MEAS_MULTIFIT_UTIL_H
#define LSST_MEAS_MULTIFIT_UTIL_H

#include <lsst/afw/image/Utils.h>
#include <ndarray.hpp>

namespace lsst {
namespace meas {
namespace multifit {


template <typename T>
static inline Eigen::Map< Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> > getMatrixView(
    ndarray::Array<T,3,3> const & array
) {
    return Eigen::Map< Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> >(
        array.getData(),
        array.template getSize<2>() * array.template getSize<1>(),
        array.template getSize<0>()
    );
}

template <typename T>
static inline Eigen::Map< Eigen::Matrix<T,Eigen::Dynamic,1> > getVectorView(
    ndarray::Array<T,2,2> const & array
) {
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

}}}

#endif //LSST_MEAS_MULTIFIT_UTIL_H
