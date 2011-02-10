#include "lsst/meas/multifit/EllipseBasis.h"
#include "lsst/meas/multifit/constants.h"

namespace lsst { namespace meas { namespace multifit {

void EllipseBasis::evaluate(
    lsst::ndarray::Array<double, 2, 1> const & matrix,
    PTR(Footprint) const & footprint,
    Ellipse const & ellipse
) const {
    detail::checkSize(
        matrix.getSize<1>(), getSize(),
        "Number of matrix columns (%d) does not match expected value (%d)."
    );

    detail::checkSize(
        matrix.getSize<0>(), footprint->getNpix(),
        "Number of matrix rows (%d) does not match expected value (%d)."
    );
    _evaluate(matrix, footprint, ellipse);
}

}}} // namespace lsst::meas::multifit
