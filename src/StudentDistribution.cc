#include "lsst/meas/multifit/StudentDistribution.h"
#include <Eigen/Cholesky>
#include <Eigen/Array>

namespace lsst { namespace meas { namespace multifit {

void StudentDistribution::draw(Random & engine, double * parameters) const {
    ensureCached();
    for (int n = 0; n < _workspace.size(); ++n) {
        _workspace[n] = engine.gaussian();
    }
    // u is drawn from chi^2_\nu distribution
    Eigen::VectorXd::Map(parameters, getSize()) 
        = getMu() 
        + ((_dof - 2.0) / engine.chisq(_dof)) * _cached->factor.part<Eigen::LowerTriangular>() * _workspace;
}

double StudentDistribution::evaluate(double const * parameters) const {
    ensureCached();
    _workspace = Eigen::VectorXd::Map(parameters, getSize()) - getMu();
    _cached->factor.part<Eigen::LowerTriangular>().solveTriangularInPlace(_workspace);
    return _cached->normalization * 
        std::pow(1.0 + _workspace.squaredNorm() / (_dof - 2.0), 0.5 * (_dof + getSize()));
}

StudentDistribution::StudentDistribution(int size, int dof) : 
    SimpleDistribution(size), _dof(dof), _workspace(getSize())
{
    if (_dof <= 2) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "StudentDistribution only supports degrees of freedom > 2."
        );
    }
}

StudentDistribution::StudentDistribution(
    Eigen::VectorXd const & mu, Eigen::MatrixXd const & sigma, int dof
) :
    SimpleDistribution(mu.size()), _dof(dof), _workspace(mu.size())
{
    if (_dof <= 2) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "StudentDistribution only supports degrees of freedom > 2."
        );
    }
    checkShape(mu, sigma);
    _mu = mu;
    _sigma = sigma;
}

StudentDistribution::StudentDistribution(StudentDistribution const & other) :
    SimpleDistribution(other), _dof(other._dof), _cached(other._cached), _workspace(other._workspace.size())
{}

StudentDistribution & StudentDistribution::operator=(StudentDistribution const & other) {
    SimpleDistribution::operator=(other);
    _dof = other._dof;
    _cached = other._cached;
    _workspace.resize(other._workspace.size());
    return *this;
}

void StudentDistribution::updateFromSamples(
    Eigen::MatrixXd const & parameters,
    Eigen::VectorXd const & weights
) {
    ensureCached();
    Eigen::MatrixXd dx = parameters.transpose() - _mu * Eigen::RowVectorXd::Ones(weights.size());
    _cached->factor.part<Eigen::LowerTriangular>().solveTriangularInPlace(dx);
    Eigen::VectorXd wgamma = weights.cwise()
        * (dx.cwise().square().colwise().sum().transpose().cwise() + (_dof - 2.0)).cwise().inverse()
        * static_cast<double>(_dof + _mu.size());
    _mu = (parameters.transpose() * wgamma).lazy();
    _mu /= wgamma.sum();
    dx = parameters.transpose() - _mu * Eigen::RowVectorXd::Ones(weights.size());
    _sigma.part<Eigen::SelfAdjoint>() = dx * wgamma.asDiagonal() * dx.transpose();
    invalidate();
}

BaseDistribution::Ptr StudentDistribution::_clone() const {
    return boost::make_shared<StudentDistribution>(*this);
}

void StudentDistribution::ensureCached() const {
    if (!_cached) {
        _cached = boost::make_shared<Cached>();
        Eigen::LLT<Eigen::MatrixXd> llt(_sigma);
        _cached->factor = llt.matrixL();
        _cached->normalization = boost::math::tgamma_ratio((_dof + getSize()) * 0.5, _dof * 0.5)
            * std::pow((_dof - 2.0) * M_PI, getSize() * 0.5)
            * std::exp(-_cached->factor.diagonal().cwise().log().sum());
    }
}

}}} // namespace lsst::meas::multifit
