#include "lsst/meas/multifit/GaussianDistribution.h"
#include <Eigen/Cholesky>
#include <Eigen/Array>

namespace lsst { namespace meas { namespace multifit {

void GaussianDistribution::draw(Random & engine, double * parameters) const {
    ensureCached();
    for (int n = 0; n < _workspace.size(); ++n) {
        _workspace[n] = engine.gaussian();
    }
    Eigen::VectorXd::Map(parameters, getSize()) = getMu() 
        + _cached->factor.part<Eigen::LowerTriangular>() * _workspace;
}

double GaussianDistribution::evaluate(double const * parameters) const {
    ensureCached();
    _workspace = Eigen::VectorXd::Map(parameters, getSize()) - getMu();
    _cached->factor.part<Eigen::LowerTriangular>().solveTriangularInPlace(_workspace);
    return _cached->normalization * std::exp(-0.5 * _workspace.squaredNorm());
}

GaussianDistribution::GaussianDistribution(int size) : 
    SimpleDistribution(size), _workspace(getSize())
{}

GaussianDistribution::GaussianDistribution(Eigen::VectorXd const & mu, Eigen::MatrixXd const & sigma) :
    SimpleDistribution(mu.size()), _workspace(mu.size())
{
    checkShape(mu, sigma);
    _mu = mu;
    _sigma = sigma;
}

GaussianDistribution::GaussianDistribution(GaussianDistribution const & other) :
    SimpleDistribution(other), _cached(other._cached), _workspace(other._workspace.size())
{}

GaussianDistribution & GaussianDistribution::operator=(GaussianDistribution const & other) {
    SimpleDistribution::operator=(other);
    _cached = other._cached;
    _workspace.resize(other._workspace.size());
    return *this;
}

void GaussianDistribution::updateFromSamples(
    Eigen::MatrixXd const & parameters,
    Eigen::VectorXd const & weights
) {
    _mu = (parameters.transpose() * weights).lazy();
    Eigen::MatrixXd dx = parameters - Eigen::VectorXd::Ones(weights.size()) * getMu().transpose();
    _sigma.part<Eigen::SelfAdjoint>() = dx.transpose() * weights.asDiagonal() * dx;
    invalidate();
}

BaseDistribution::Ptr GaussianDistribution::_clone() const {
    return boost::make_shared<GaussianDistribution>(*this);
}

void GaussianDistribution::ensureCached() const {
    if (!_cached) {
        _cached = boost::make_shared<Cached>();
        Eigen::LLT<Eigen::MatrixXd> llt(_sigma);
        _cached->factor = llt.matrixL();
        _cached->normalization = std::exp(-_cached->factor.diagonal().cwise().log().sum())
            * std::pow(2.0 * M_PI, -0.5 * getSize());
    }
}

}}} // namespace lsst::meas::multifit
