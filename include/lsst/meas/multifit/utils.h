#ifndef LSST_MEAS_MULTIFIT_UTILS_H
#define LSST_MEAS_MULTIFIT_UTILS_H

#include "lsst/afw/image/Utils.h"
#include <ndarray.hpp>

namespace lsst {
namespace meas {
namespace multifit {

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

#endif //LSST_MEAS_MULTIFIT_UTILS_H
