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

lsst::ndarray::Array<Pixel const,1,1> ModelBasis::getIntegration() const {
    if (_multipoleMatrix.getData() == 0) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "getIntegration() not implemented for this basis."
        );
    }
    return _multipoleMatrix[0];
}

MultipoleMatrix ModelBasis::getMultipoleMatrix() const {
    if (_multipoleMatrix.getData() == 0 || _multipoleMatrix.getSize<0>() != 6) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "getMultipoleMatrix() not implemented for this basis."
        );
    }
    return MultipoleMatrix(_multipoleMatrix);
}

ModelBasis::Ptr ModelBasis::convolve(CONST_PTR(LocalPsf) const & psf) const {
    throw LSST_EXCEPT(lsst::pex::exceptions::LogicErrorException,
                      "ModelBasis subclass does not support convolution.");
}

void ModelBasis::attachConstraint(
    lsst::ndarray::Array<Pixel const,2,1> const & matrix,
    lsst::ndarray::Array<Pixel const,1,1> const & vector
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

void ModelBasis::attachMultipoleMatrix(
    lsst::ndarray::Array<Pixel const,2,2> const & matrix
) {
    detail::checkSize(
        matrix.getSize<0>(), 6,
        "Number of rows of multipole matrix (%d) must be %d."
    );
    detail::checkSize(
        matrix.getSize<1>(), _size,
        "Number of columns of multipole matrix (%d) must match basis size (%d)."
    );
    _multipoleMatrix = matrix;
}

}}} // namespace lsst::meas::multifit
