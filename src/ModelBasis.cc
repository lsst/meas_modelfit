#include "lsst/meas/multifit/ModelBasis.h"
#include "lsst/meas/multifit/constants.h"

namespace lsst { namespace meas { namespace multifit {

void ModelBasis::evaluate(
    lsst::ndarray::Array<Pixel, 2, 1> const & matrix,
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

void ModelBasis::evaluateRadialProfile(
    lsst::ndarray::Array<Pixel,2,1> const & profile,
    lsst::ndarray::Array<Pixel const,1,1> const & radii
) const {
    detail::checkSize(
        profile.getSize<1>(), getSize(),
        "Number of matrix columns (%d) does not match expected value (%d)."
    );

    detail::checkSize(
        profile.getSize<0>(), radii.getSize<0>(),
        "Number of matrix rows (%d) does not match expected value (%d)."
    );
    _evaluateRadialProfile(profile, radii);
}

void ModelBasis::integrate(lsst::ndarray::Array<Pixel, 1, 1> const & vector) const {
    detail::checkSize(
        vector.getSize<0>(), getSize(),
        "Number of vector elements (%d) does not match expected value (%d)."
    );
    _integrate(vector);
}

void ModelBasis::evaluateMultipoleMatrix(lsst::ndarray::Array<Pixel, 2, 1> const & matrix) const {
    detail::checkSize(
        matrix.getSize<1>(), getSize(),
        "Number of coefficients in matrix (%d) does not match expected value (%d)."
    );
    detail::checkSize(
        matrix.getSize<0>(), 6,
        "Incorrect number of rows (%d) in multipole matrix (expected %d)."
    );
    _evaluateMultipoleMatrix(matrix);
}

ModelBasis::Ptr ModelBasis::convolve(CONST_PTR(LocalPsf) const & psf) const {
    throw LSST_EXCEPT(lsst::pex::exceptions::LogicErrorException,
                      "ModelBasis subclass does not support convolution.");
}

void ModelBasis::attachConstraint(
    lsst::ndarray::Array<Pixel,2,1> const & matrix,
    lsst::ndarray::Array<Pixel,1,1> const & vector
) {
    detail::checkSize(
        matrix.getSize<0>(), vector.getSize<0>(),
        "Number of constraints in matrix (%d) do not match number of constraints in vector (%d)."
    );
    detail::checkSize(
        matrix.getSize<1>(), _size,
        "Incorrect number of columns (%d) in constraint matrix (expected %d)."
    );
    _constraintMatrix = matrix;
    _constraintVector = vector;
}

}}} // namespace lsst::meas::multifit
