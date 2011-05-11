#include "lsst/meas/multifit/ModelBasis.h"
#include "lsst/meas/multifit/constants.h"

namespace lsst { namespace meas { namespace multifit {

void ModelBasis::evaluate(
    lsst::ndarray::Array<double, 2, 1> const & matrix,
    CONST_PTR(Footprint) const & footprint,
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

void ModelBasis::integrate(lsst::ndarray::Array<double, 1, 1> const & vector) const {
    detail::checkSize(
        vector.getSize<0>(), getSize(),
        "Number of vector elements (%d) does not match expected value (%d)."
    );
    _integrate(vector);
}

ModelBasis::Ptr ModelBasis::convolve(CONST_PTR(LocalPsf) const & psf) const {
    throw LSST_EXCEPT(lsst::pex::exceptions::LogicErrorException,
                      "ModelBasis subclass does not support convolution.");
}

}}} // namespace lsst::meas::multifit
