#include "lsst/meas/multifit/constants.h"
#include "lsst/meas/multifit/BaseEvaluator.h"
#include "lsst/ndarray/eigen.h"
#include <limits>

namespace lsst { namespace meas { namespace multifit {

void BaseEvaluator::evaluateModelMatrix(
    ndarray::Array<Pixel,2,2> const & matrix,
    ndarray::Array<double const,1,1> const & param
) const {
    detail::checkSize(
        matrix.getSize<0>(), getPixelCount(),
        "Number of matrix rows (%d) does not match expected value (%d)."
    );
    detail::checkSize(
        matrix.getSize<1>(), getCoefficientCount(),
        "Number of matrix columns (%d) does not match expected value (%d)."
    );
    detail::checkSize(
        param.getSize<0>(), getParameterCount(),
        "Parameter vector size (%d) does not match expected value (%d)."
    );
    _evaluateModelMatrix(matrix, param);
}

void BaseEvaluator::evaluateModelMatrixDerivative(
    ndarray::Array<Pixel,3,3> const & derivative,
    ndarray::Array<double const,1,1> const & param
) const {
    detail::checkSize(
        derivative.getSize<0>(), getParameterCount(),
        "Size of derivative array first dimension (%d) does not match expected value (%d)."
    );
    detail::checkSize(
        derivative.getSize<1>(), getPixelCount(),
        "Size of derivative array second dimension (%d) does not match expected value (%d)."
    );
    detail::checkSize(
        derivative.getSize<2>(), getCoefficientCount(),
        "Size of derivative array third (%d) does not match expected value (%d)."
    );
    detail::checkSize(
        param.getSize<0>(), getParameterCount(),
        "Parameter vector size (%d) does not match expected value (%d)."
    );
    ndarray::Array<Pixel,2,2> modelMatrix = ndarray::allocate(getPixelCount(), getCoefficientCount());
    _evaluateModelMatrix(modelMatrix, param);
    _evaluateModelMatrixDerivative(derivative, modelMatrix, param);
}

void BaseEvaluator::_evaluateModelMatrixDerivative(
    ndarray::Array<Pixel,3,3> const & derivative,
    ndarray::Array<Pixel const,2,2> const & fiducial,
    ndarray::Array<double const,1,1> const & param
) const {
    static double const epsilon = std::sqrt(
        std::numeric_limits<double>::epsilon()
    );
    ndarray::Array<double,1,1> parameters(ndarray::copy(param));
    for (int n = 0; n < getParameterCount(); ++n) {
        parameters[n] += epsilon;
        _evaluateModelMatrix(derivative[n], parameters);
        derivative[n] -= fiducial;
        derivative[n] /= epsilon;
        parameters[n] -= epsilon;
    }
}

void BaseEvaluator::writeInitialParameters(
    ndarray::Array<double,1,1> const & param
) const {
    detail::checkSize(
        param.getSize<0>(), getParameterCount(),
        "Parameter vector size (%d) does not match expected value (%d)."
    );
    _writeInitialParameters(param);
}

}}} // namespace lsst::meas::multifit
