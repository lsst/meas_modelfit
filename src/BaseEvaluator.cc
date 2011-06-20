#include "lsst/meas/multifit/BaseEvaluator.h"
#include "lsst/ndarray/eigen.h"
#include <Eigen/Cholesky>

namespace lsst { namespace meas { namespace multifit {

CoefficientPrior::ConstPtr BaseEvaluator::evaluate(
    ndarray::Array<Pixel,2,2> const & matrix,
    ndarray::Array<double const,1,1> const & parameters
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
        parameters.getSize<0>(), getParameterCount(),
        "Parameter vector size (%d) does not match expected value (%d)."
    );
    return _evaluate(matrix, parameters);
}

void BaseEvaluator::writeInitialParameters(
    ndarray::Array<double,1,1> const & parameters
) const {
    detail::checkSize(
        parameters.getSize<0>(), getParameterCount(),
        "Parameter vector size (%d) does not match expected value (%d)."
    );
    _writeInitialParameters(parameters);
}

double BaseEvaluator::clipToBounds(ndarray::Array<double,1,1> const & parameters) const {
    detail::checkSize(
        parameters.getSize<0>(), getParameterCount(),
        "Parameter vector size (%d) does nto match expected value (%d)."
    );
    return _clipToBounds(parameters);
}

bool BaseEvaluator::checkBounds(ndarray::Array<double const,1,1> const & parameters) const {
    detail::checkSize(
        parameters.getSize<0>(), getParameterCount(),
        "Parameter vector size (%d) does nto match expected value (%d)."
    );
    return _checkBounds(parameters);
}

double BaseEvaluator::integrate(
    Random & engine,
    ndarray::Array<Pixel,2,2> const & coefficients,
    ndarray::Array<Pixel,1,1> const & weights,
    ndarray::Array<double const,1,1> const & parameters
) const {
    ndarray::Array<Pixel,2,2> modelMatrix(ndarray::allocate(getPixelCount(), getCoefficientCount()));
    CoefficientPrior::ConstPtr prior = _evaluate(modelMatrix, parameters);
    double result = prior->integrate(engine, coefficients, weights, modelMatrix, getDataVector());
    result *= std::exp(getLogPixelErrorSum());
    return result;
}

}}} // namespace lsst::meas::multifit
