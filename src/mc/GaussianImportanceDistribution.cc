#include "lsst/meas/multifit/mc/GaussianImportanceDistribution.h"
#include <Eigen/Cholesky>
#include <Eigen/Array>
#include "lsst/ndarray/eigen.h"

namespace lsst { namespace meas { namespace multifit { namespace mc {

void GaussianImportanceDistribution::draw(
    Random & engine,
    ndarray::Array<double,2,2> const & parameters,
    ndarray::Array<double,1,1> const & importance
) const {
    for (ndarray::Array<double,2,2>::Iterator i = parameters.begin(); i != parameters.end(); ++i) {
        for (int n = 0; n < parameters.getSize<0>(); ++n) {
            _workspace[n] = engine.gaussian();
        }
        ndarray::viewAsEigen(*i) = _factor.part<Eigen::LowerTriangular>() * _workspace;
        ndarray::viewAsEigen(*i) += _mean;
    }
    importance.deep() = 0.0;
    evaluate(parameters, importance);
}

void GaussianImportanceDistribution::evaluate(
    ndarray::Array<double const,2,2> const & parameters,
    ndarray::Array<double,1,1> const & output,
    double factor
) const {
    ndarray::Array<double const,2,2>::Iterator pi = parameters.begin();
    ndarray::Array<double,1,1>::Iterator oi = output.begin();
    for (; pi != parameters.end(); ++pi, ++oi) {
        _workspace = ndarray::viewAsEigen(*pi) - _mean;
        _factor.part<Eigen::LowerTriangular>().solveTriangularInPlace(_workspace);
        *oi += factor * _normalization * std::exp(-0.5 * _workspace.squaredNorm());
    }
}

ImportanceDistribution::Ptr GaussianImportanceDistribution::adapt(
    ndarray::Array<double const,2,1> const & parameters,
    ndarray::Array<double const,1,1> const & weights
) const {
    ndarray::EigenView<double const,2,1> p(parameters);
    ndarray::EigenView<double const,1,1> w(weights);
    Eigen::VectorXd mean = p.transpose() * w;
    Eigen::MatrixXd dx = p - Eigen::VectorXd::Ones(w.size()) * _mean.transpose();
    Eigen::MatrixXd sigma(mean.size(), mean.size());
    sigma.part<Eigen::SelfAdjoint>() = dx.transpose() * w.asDiagonal() * dx;
    return ImportanceDistribution::Ptr(new GaussianImportanceDistribution(mean, sigma));
}

GaussianImportanceDistribution::Ptr GaussianImportanceDistribution::make(
    ndarray::Array<double const,1,1> const & mean,
    ndarray::Array<double const,2,2> const & sigma
) {
    detail::checkSize(
        sigma.getSize<0>(), mean.getSize<0>(),
        "Number of rows of covariance matrix (%d) does not match size of mean vector (%d)."
    );
    detail::checkSize(
        sigma.getSize<1>(), mean.getSize<0>(),
        "Number of columns of covariance matrix (%d) does not match size of mean vector (%d)."
    );
    Eigen::VectorXd m(ndarray::viewAsEigen(mean));
    Eigen::MatrixXd s(ndarray::viewAsEigen(sigma));
    return GaussianImportanceDistribution::Ptr(new GaussianImportanceDistribution(m, s));
}

GaussianImportanceDistribution::Ptr GaussianImportanceDistribution::make(
    ndarray::Array<double const,1,1> const & mean,
    ndarray::Array<double const,1,1> const & sigmaDiagonal
) {
    detail::checkSize(
        sigmaDiagonal.getSize<0>(), mean.getSize<0>(),
        "Number of elements in covariance matrix diagonal (%d) does not match size of mean vector (%d)."
    );
    Eigen::VectorXd m(ndarray::viewAsEigen(mean));
    Eigen::MatrixXd s = Eigen::MatrixXd::Zero(m.size(), m.size());
    s.diagonal() = ndarray::viewAsEigen(sigmaDiagonal);
    return GaussianImportanceDistribution::Ptr(new GaussianImportanceDistribution(m, s));
}

GaussianImportanceDistribution::GaussianImportanceDistribution(
    Eigen::VectorXd & mean, Eigen::MatrixXd & sigma
) : _normalization(0.0), _mean(), _factor(), _workspace(mean.size())
{
    _mean.swap(mean);
    Eigen::LLT<Eigen::MatrixXd> cholesky(sigma);
    _factor.swap(sigma);
    _factor = cholesky.matrixL();
    _normalization = std::exp(-_factor.diagonal().cwise().log().sum()) 
        * std::pow(2.0 * M_PI, -0.5 * _factor.rows());
}

}}}} // namespace lsst::meas::multifit::mc
