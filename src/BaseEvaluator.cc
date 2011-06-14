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

double BaseEvaluator::integrate(
    Random & engine,
    ndarray::Array<Pixel,2,2> const & coefficients,
    ndarray::Array<Pixel,1,1> const & weights,
    ndarray::Array<double const,1,1> const & parameters
) const {
    detail::checkSize(
        parameters.getSize<0>(), getParameterCount(),
        "Parameter vector size (%d) does nto match expected value (%d)."
    );
    detail::checkSize(
        coefficients.getSize<1>(), getCoefficientCount(),
        "Number of coefficient array columns (%d) does not match expected value (%d)."
    );
    detail::checkSize(
        coefficients.getSize<0>(), weights.getSize<0>(),
        "Number of coefficient array rows (%d) does not match weights array size (%d)."
    );
    ndarray::Array<Pixel,2,2> modelMatrix(ndarray::allocate(getPixelCount(), getCoefficientCount()));
    CoefficientPrior::ConstPtr prior = _evaluate(modelMatrix, parameters);
    Eigen::MatrixXd h(getCoefficientCount(), getCoefficientCount());
    h.part<Eigen::SelfAdjoint>() =
        ndarray::viewAsTransposedEigen(modelMatrix) * ndarray::viewAsEigen(modelMatrix);
    h += prior->getGaussianMatrix();
    Eigen::VectorXd mu = ndarray::viewAsTransposedEigen(modelMatrix) * ndarray::viewAsEigen(getDataVector());
    mu += prior->getGaussianMatrix() * prior->getGaussianVector();
    Eigen::LLT<Eigen::MatrixXd> cholesky(h);
    if (!cholesky.isPositiveDefinite()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::RuntimeErrorException,
            "Cholesky decomposition failure in Monte Carlo integration."
        );
    }
    cholesky.solveInPlace(mu);
    Eigen::VectorXd r = ndarray::viewAsEigen(modelMatrix) * mu - ndarray::viewAsEigen(getDataVector());
    long double s = 0.5 * r.squaredNorm();
    r = mu - prior->getGaussianVector();
    s += 0.5 * r.dot(prior->getGaussianMatrix() * r);
    s += getLogPixelErrorSum();
    s -= (cholesky.matrixL().diagonal() / (2.0 * M_PI)).cwise().log().sum();
    long double result = 0.0;
    ndarray::Array<Pixel,2,2>::Iterator ci;
    ndarray::Array<Pixel,1,1>::iterator wi;
    for (ci = coefficients.begin(), wi = weights.begin(); wi != weights.end(); ++ci, ++wi) {
        for (ndarray::Array<Pixel,2,2>::Reference::Iterator cj = ci->begin(); cj != ci->end(); ++cj) {
            *cj = engine.gaussian();
        }
        ndarray::EigenView<double,1,1> coeff(*ci);
        cholesky.matrixL().solveTriangularInPlace(coeff);
        coeff += mu;
        result += (*wi) = (*prior)(*ci);
    }
    ndarray::viewAsEigen(weights) /= result;
    result *= std::exp(-s) / weights.getSize<0>();
    return result;
}

}}} // namespace lsst::meas::multifit
