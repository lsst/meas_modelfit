#include "lsst/meas/multifit/GaussianDistribution.h"
#include <Eigen/Cholesky>
#include <Eigen/Array>
#include "lsst/ndarray/eigen.h"

namespace lsst { namespace meas { namespace multifit {

namespace {

double computeNormalization(Eigen::MatrixXd const & factor) {
    return std::exp(-factor.diagonal().cwise().log().sum()) * std::pow(2.0 * M_PI, -0.5 * factor.rows());
}

} // anonymous


void GaussianDistribution::draw(Random & engine, double * parameters) const {
    ensureCached();
    for (int n = 0; n < _workspace.size(); ++n) {
        _workspace[n] = engine.gaussian();
    }
    Eigen::VectorXd::Map(parameters, getDimensionality()) = getMu() 
        + _cached->factor.part<Eigen::LowerTriangular>() * _workspace;
}

double GaussianDistribution::evaluate(double const * parameters) const {
    ensureCached();
    _workspace = Eigen::VectorXd::Map(parameters, getDimensionality()) - getMu();
    _cached->factor.part<Eigen::LowerTriangular>().solveTriangularInPlace(_workspace);
    return _cached->normalization * std::exp(-0.5 * _workspace.squaredNorm());
}

int GaussianDistribution::getNestedDimensionality() const {
    return (_nested) ? _nested->getDimensionality() : 0;
}

void GaussianDistribution::convertUnifiedToNested(int nx) {
    int ny = getDimensionality() - nx;
    if (nx < 0 || ny < 0) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            (boost::format("Cannot create nested distribution with negative dimensions (%d, %d).")
             % nx % ny).str()
        );
    }
    ensureCached();
    Eigen::MatrixXd fisher = Eigen::MatrixXd::Identity(nx + ny, nx + ny);
    _cached->factor.part<Eigen::LowerTriangular>().solveTriangularInPlace(fisher);
    _cached->factor.part<Eigen::LowerTriangular>().transpose().solveTriangularInPlace(fisher);
    _nested.reset(new GaussianDistribution(ny));
    _nested->_mu = _mu.segment(nx, ny);
    Eigen::LLT<Eigen::MatrixXd> xCholesky(fisher.block(0, 0, nx, nx));
    _sigma = Eigen::MatrixXd::Identity(nx, nx);
    _cached->factor = Eigen::MatrixXd::Identity(nx, nx);
    xCholesky.matrixL().solveTriangularInPlace(_cached->factor);
    _cached->normalization = computeNormalization(_cached->factor);
    xCholesky.solveInPlace(_sigma);
    Eigen::LLT<Eigen::MatrixXd> yCholesky(fisher.block(nx, nx, ny, ny));
    _nested->_cached = boost::make_shared<Cached>();
    _nested->_cached->factor = Eigen::MatrixXd::Identity(ny, ny);
    yCholesky.matrixL().solveTriangularInPlace(_nested->_cached->factor);
    _nested->_cached->normalization = computeNormalization(_nested->_cached->factor);
    yCholesky.solveInPlace(_nested->_sigma);
    double threshold = fisher.maxCoeff() * std::sqrt(std::numeric_limits<double>::epsilon());
    if (!fisher.block(0, nx, nx, ny).isZero(threshold)) {
        _nestedConditional.reset(new Eigen::MatrixXd(fisher.block(0, nx, nx, ny)));
    }
    _mu = _mu.segment(0, nx);
    _dimensionality = nx;
}

void GaussianDistribution::convertNestedToUnified() {
    if (!_nested) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Distribution is not currently nested."
        );
    }
    int nx = getDimensionality();
    int ny = _nested->getDimensionality();
    ensureCached();
    Eigen::MatrixXd fisher = Eigen::MatrixXd::Identity(nx + ny, nx + ny);
    _cached->factor.part<Eigen::LowerTriangular>()
        .solveTriangularInPlace(fisher.block(0, 0, nx, nx));
    _cached->factor.part<Eigen::LowerTriangular>()
        .transpose().solveTriangularInPlace(fisher.block(0, 0, nx, nx));
    _nested->ensureCached();
    _nested->_cached->factor.part<Eigen::LowerTriangular>()
        .solveTriangularInPlace(fisher.block(nx, nx, ny, ny));
    _nested->_cached->factor.part<Eigen::LowerTriangular>()
        .transpose().solveTriangularInPlace(fisher.block(nx, nx, ny, ny));
    if (_nestedConditional) {
        fisher.block(0, nx, nx, ny) = *_nestedConditional;
        fisher.block(nx, 0, ny, nx) = _nestedConditional->transpose();
        _nestedConditional.reset();
    }
    Eigen::LLT<Eigen::MatrixXd> cholesky(fisher);
    _cached->factor = Eigen::MatrixXd::Identity(nx + ny, nx + ny);
    cholesky.matrixL().solveTriangularInPlace(_cached->factor);
    _cached->normalization = computeNormalization(_cached->factor);
    _sigma = Eigen::MatrixXd::Identity(nx + ny, nx + ny);
    cholesky.solveInPlace(_sigma);
    Eigen::VectorXd mu(nx + ny);
    mu.segment(0, nx) = _mu;
    mu.segment(nx, ny) = _nested->_mu;
    _mu.swap(mu);
    _nested.reset();
    _dimensionality = nx + ny;
}

GaussianDistribution::GaussianDistribution(int dimensionality) : 
    SimpleDistribution(dimensionality), _workspace(getDimensionality())
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
    ndarray::Array<double const,2,1> const & parameters,
    ndarray::Array<double const,1,1> const & weights
) {
    ndarray::EigenView<double const,2,1> p(parameters);
    ndarray::EigenView<double const,1,1> w(weights);
    _mu = (p.transpose() * w).lazy();
    Eigen::MatrixXd dx = p - Eigen::VectorXd::Ones(w.size()) * getMu().transpose();
    _sigma.part<Eigen::SelfAdjoint>() = dx.transpose() * w.asDiagonal() * dx;
    invalidate();
}

BaseDistribution::Ptr GaussianDistribution::_clone() const {
    return boost::make_shared<GaussianDistribution>(*this);
}

BaseDistribution::Ptr GaussianDistribution::_evaluateNested(double const * parameters) const {
    GaussianDistribution::Ptr result;
    if (_nested) return result;
    result = _nested->clone();
    _updateNested(*result, parameters);
    return result;   
}

void GaussianDistribution::_updateNested(BaseDistribution & nested, double const * parameters) const {
    if (_nestedConditional) {
        GaussianDistribution & result = static_cast<GaussianDistribution &>(nested);
        _workspace = Eigen::VectorXd::Map(parameters, getDimensionality()) - _mu;
        _nested->_workspace = (_nestedConditional->transpose() * _workspace).lazy();
        result._workspace = (result._sigma * _nested->_workspace).lazy();
        result._mu = _nested->_mu - result._workspace;
    }
}

void GaussianDistribution::ensureCached() const {
    if (!_cached) {
        _cached = boost::make_shared<Cached>();
        Eigen::LLT<Eigen::MatrixXd> llt(_sigma);
        _cached->factor = llt.matrixL();
        _cached->normalization = computeNormalization(_cached->factor);
    }
}

}}} // namespace lsst::meas::multifit
