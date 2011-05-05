#include "lsst/meas/multifit/constants.h"
#include "lsst/meas/multifit/BaseEvaluator.h"
#include "lsst/ndarray/eigen.h"
#include <limits>

namespace lsst { namespace meas { namespace multifit {

bool BaseEvaluator::evaluateModelMatrix(
    ndarray::Array<double,2,2> const & matrix,
    ndarray::Array<double const,1,1> const & param
) const {
    detail::checkSize(
        matrix.getSize<0>(), getDataSize(),
        "Number of matrix rows (%d) does not match expected value (%d)."
    );
    detail::checkSize(
        matrix.getSize<1>(), getCoefficientSize(),
        "Number of matrix columns (%d) does not match expected value (%d)."
    );
    detail::checkSize(
        param.getSize<0>(), getParameterSize(),
        "Parameter vector size (%d) does not match expected value (%d)."
    );
    return _evaluateModelMatrix(matrix, param);
}

bool BaseEvaluator::evaluateModelMatrixDerivative(
    ndarray::Array<double,3,3> const & derivative,
    ndarray::Array<double const,1,1> const & param
) const {
    detail::checkSize(
        derivative.getSize<0>(), getParameterSize(),
        "Size of derivative array first dimension (%d) does not match expected value (%d)."
    );
    detail::checkSize(
        derivative.getSize<1>(), getDataSize(),
        "Size of derivative array second dimension (%d) does not match expected value (%d)."
    );
    detail::checkSize(
        derivative.getSize<2>(), getCoefficientSize(),
        "Size of derivative array third (%d) does not match expected value (%d)."
    );
    detail::checkSize(
        param.getSize<0>(), getParameterSize(),
        "Parameter vector size (%d) does not match expected value (%d)."
    );
    ndarray::Array<double,2,2> modelMatrix = ndarray::allocate(getDataSize(), getCoefficientSize());
    _evaluateModelMatrix(modelMatrix, param);
    return _evaluateModelMatrixDerivative(derivative, modelMatrix, param);
}

bool BaseEvaluator::_evaluateModelMatrixDerivative(
    ndarray::Array<double,3,3> const & derivative,
    ndarray::Array<double const,2,2> const & fiducial,
    ndarray::Array<double const,1,1> const & param
) const {
    static double const epsilon = std::sqrt(
        std::numeric_limits<double>::epsilon()
    );
    ndarray::Array<double, 1, 1> parameters(ndarray::copy(param));
    for (int n = 0; n < _parameterSize; ++n) {
        parameters[n] += epsilon;
        if (!_evaluateModelMatrix(derivative[n], parameters)) return false;
        derivative[n] -= fiducial;
        derivative[n] /= epsilon;
        parameters[n] -= epsilon;
    }
    return true;
}

void BaseEvaluator::writeInitialParameters(
    ndarray::Array<double,1,1> const & param
) const {
    detail::checkSize(
        param.getSize<0>(), getParameterSize(),
        "Parameter vector size (%d) does not match expected value (%d)."
    );
    _writeInitialParameters(param);
}

}}} // namespace lsst::meas::multifit
