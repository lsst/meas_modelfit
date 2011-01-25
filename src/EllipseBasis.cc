#include "lsst/meas/multifit/EllipseBasis.h"

namespace lsst { namespace meas { namespace multifit {

void EllipseBasis::evaluate(
    FootprintMatrix const & matrix,
    lsst::afw::geom::Ellipse const & ellipse
) const {
    checkSizes(
        matrix.getArray().getSize<1>(), getSize(),
        "Number of matrix columns (%d) does not match expected value (%d)."
    );
    _evaluate(matrix, ellipse);
}

}}} // namespace lsst::meas::multifit
