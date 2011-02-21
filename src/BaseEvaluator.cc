#include "lsst/meas/multifit/constants.h"
#include "lsst/meas/multifit/BaseEvaluator.h"
#include "lsst/ndarray/eigen.h"
#include <limits>

namespace lsst { namespace meas { namespace multifit {

void BaseEvaluator::evaluateModelMatrix(
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
    _evaluateModelMatrix(matrix, param);
}

void BaseEvaluator::evaluateModelDerivative(
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
    _evaluateModelDerivative(derivative, param);
}

void BaseEvaluator::_evaluateModelDerivative(
    ndarray::Array<double,3,3> const & derivative,
    ndarray::Array<double const,1,1> const & param
) const {
    static double const epsilon = std::sqrt(
        std::numeric_limits<double>::epsilon()
    );

    ndarray::Array<double,2,2> fiducial(
        ndarray::allocate(
            ndarray::makeVector(getDataSize(), getCoefficientSize())
        )
    );
    _evaluateModelMatrix(fiducial, param);
    ndarray::Array<double, 1, 1> parameters(ndarray::copy(param));
    for (int n = 0; n < _parameterSize; ++n) {
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
        param.getSize<0>(), getParameterSize(),
        "Parameter vector size (%d) does not match expected value (%d)."
    );
    _writeInitialParameters(param);
}

}}} // namespace lsst::meas::multifit
