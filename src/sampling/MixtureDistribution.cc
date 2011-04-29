#include "lsst/meas/multifit/sampling/MixtureDistribution.h"
#include "lsst/ndarray/eigen.h"
#include "lsst/pex/exceptions.h"

#include <boost/math/special_functions/gamma.hpp>
#include <Eigen/Array>
#include <Eigen/Cholesky>

namespace lsst { namespace meas { namespace multifit { namespace sampling {

namespace {

double evaluateGaussian(
    Eigen::VectorXd const & mu, Eigen::MatrixXd const & sigma,
    ndarray::Array<double const,1,1> const & parameters, 
    Eigen::VectorXd & workspace
) {
    workspace = ndarray::viewAsEigen(parameters) - mu;
    sigma.part<Eigen::LowerTriangular>().solveTriangularInPlace(workspace);
    return std::exp(
        - 0.5 * workspace.squaredNorm() 
        - sigma.diagonal().cwise().log().sum()
    );
}

void drawGaussian(
    Random & engine,
    Eigen::VectorXd const & mu, Eigen::MatrixXd const & sigma,
    ndarray::Array<double,1,1> const & parameters,
    Eigen::VectorXd & workspace
) {
    for (int n = 0; n < workspace.size(); ++n) {
        workspace[n] = engine.gaussian();
    }
    ndarray::viewAsEigen(parameters) = mu + sigma.part<Eigen::LowerTriangular>() * workspace;
} 

void updateGaussian(
    ConstTable const & table, double constant, 
    double & normalization, Eigen::VectorXd & mu, Eigen::MatrixXd & sigma
) {
    Eigen::VectorXd rho(table.getTableSize());
    Eigen::VectorXd tau(table.getTableSize());
    Eigen::VectorXd workspace(table.getParameterSize());
    double alpha = 0.0;
    for (int i = 0; i < table.getTableSize(); ++i) {
        ConstRecord record = table[i];
        rho[i] = normalization * constant
            * evaluateGaussian(mu, sigma, record.parameters, workspace)
            / record.importance;
        tau[i] = record.weight * rho[i] / table.getWeightSum();
        alpha += tau[i];
    }
    tau /= alpha;
    mu.setZero();
    for (int i = 0; i < table.getTableSize(); ++i) {
        ConstRecord record = table[i];
        mu += tau[i] * ndarray::viewAsEigen(record.parameters);
    }
    sigma.setZero();
    for (int i = 0; i < table.getTableSize(); ++i) {
        ConstRecord record = table[i];
        workspace = ndarray::viewAsEigen(record.parameters) - mu;
        sigma += tau[i] * workspace * workspace.transpose();
    }
    normalization = alpha;
    Eigen::LLT<Eigen::MatrixXd> cholesky(sigma);
    sigma = cholesky.matrixL();
}

double evaluateStudent(
    int dof, Eigen::VectorXd const & mu, Eigen::MatrixXd const & sigma,
    ndarray::Array<double const,1,1> const & parameters, 
    Eigen::VectorXd & workspace
) {
    workspace = ndarray::viewAsEigen(parameters) - mu;
    sigma.part<Eigen::LowerTriangular>().solveTriangularInPlace(workspace);
    return std::exp(-sigma.diagonal().cwise().log().sum())
        * std::pow(1.0 + workspace.squaredNorm() / dof, -0.5*(dof + workspace.size()));
}

void drawStudent(
    Random & engine,
    int dof, Eigen::VectorXd const & mu, Eigen::MatrixXd const & sigma,
    ndarray::Array<double,1,1> const & parameters,
    Eigen::VectorXd & workspace
) {
    for (int n = 0; n < workspace.size(); ++n) {
        workspace[n] = engine.gaussian();
    }
    ndarray::viewAsEigen(parameters) = mu 
        + (sigma.part<Eigen::LowerTriangular>() * workspace)
        * std::sqrt(0.5 * dof / boost::math::gamma_p_inv(0.5 * dof, engine.uniform()));
}

void updateStudent(
    ConstTable const & table, int dof, double constant,
    double & normalization, Eigen::VectorXd & mu, Eigen::MatrixXd & sigma
) {
    Eigen::VectorXd rho(table.getTableSize());
    Eigen::VectorXd gamma(table.getTableSize());
    Eigen::VectorXd tau(table.getTableSize());
    Eigen::VectorXd workspace(table.getParameterSize());
    double alpha = 0.0;
    double beta = 0.0;
    for (int i = 0; i < table.getTableSize(); ++i) {
        ConstRecord record = table[i];
        rho[i] = normalization * constant
            * evaluateStudent(dof, mu, sigma, record.parameters, workspace) / record.importance;
        gamma[i] = (dof + mu.size()) / (dof + workspace.squaredNorm());
        tau[i] = record.weight * rho[i] / table.getWeightSum();
        alpha += tau[i];
        beta += tau[i] * gamma[i];
    }
    normalization = alpha;
    mu.setZero();
    for (int i = 0; i < table.getTableSize(); ++i) {
        ConstRecord record = table[i];
        mu += tau[i] * gamma[i] * ndarray::viewAsEigen(record.parameters) / beta;
    }
    sigma.setZero();
    for (int i = 0; i < table.getTableSize(); ++i) {
        ConstRecord record = table[i];
        workspace = ndarray::viewAsEigen(record.parameters) - mu;
        sigma += tau[i] * gamma[i] * workspace * workspace.transpose() / alpha;
    }
    Eigen::LLT<Eigen::MatrixXd> cholesky(sigma);
    sigma = cholesky.matrixL();
}

} // anonymous

MixtureDistribution::MixtureDistribution(ComponentList const & components, int dof) :
    _dof(dof), _constant(0.0), _components(components)
{
    int const p = _components.back().getMu().size();
    if (_dof < 0) {
        _constant = std::pow(2.0 * M_PI, -p * 0.5);
    } else {
        _constant = boost::math::tgamma_ratio((_dof + p)*0.5, _dof*0.5)
            * std::pow(_dof * M_PI, p * 0.5);
    }
    double total = 0.0;
    for (ComponentList::iterator j = _components.begin(); j != _components.end(); ++j) {
        total += j->getNormalization();
    }
    for (ComponentList::iterator j = _components.begin(); j != _components.end(); ++j) {
        j->_normalization /= total;
    }
}

void MixtureDistribution::draw(Table const & table, Random & engine) const {
    if (table.getParameterSize() != getParameterSize()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LengthErrorException,
            "Number of parameters in table does not match number of parameters in distribution"
        );
    }
    Eigen::VectorXd workspace(getParameterSize());
    table.clear();
    for (int i = 0; i < table.getTableSize(); ++i) {
        Record record = table[i];
        double u = engine.uniform();
        double v = 0.0;
        ComponentList::const_iterator j = _components.begin();
        for (; j != _components.end(); ++j) {
            v += j->getNormalization();
            if (v >= u) break;
        }
        if (j == _components.end()) --j; // only needed because of round-off error
        if (_dof < 0) {
            drawGaussian(engine, j->getMu(), j->getSigma(), record.parameters, workspace);
        } else {
            drawStudent(engine, _dof, j->getMu(), j->getSigma(), record.parameters, workspace);
        }
        // Now evaluate the probability at that point.
        for (j = _components.begin(); j != _components.end(); ++j) {
            double q;
            if (_dof < 0) {
                q = evaluateGaussian(j->getMu(), j->getSigma(), record.parameters, workspace);
            } else {
                q = evaluateStudent(_dof, j->getMu(), j->getSigma(), record.parameters, workspace);
            }
            record.importance += j->getNormalization() * _constant * q;
        }
    }
}

void MixtureDistribution::update(ConstTable const & table) {
    for (ComponentList::iterator j = _components.begin(); j != _components.end(); ++j) {
        if (_dof < 0) {
            updateGaussian(table, _constant, j->_normalization, j->_mu, j->_sigma);
        } else {
            updateStudent(table, _dof, _constant, j->_normalization, j->_mu, j->_sigma);
        }
    }
}

MixtureDistribution MixtureDistribution::scatter(
    Random & engine,
    int nComponents, double fraction, int dof,
    lsst::ndarray::Array<double const,1,1> const & mean,
    lsst::ndarray::Array<double const,2,2> const & covariance
) {
    MixtureDistribution::ComponentList components;
    Eigen::LLT<Eigen::MatrixXd> cholesky(ndarray::viewAsEigen(covariance));
    Eigen::VectorXd workspace(mean.getSize<0>());
    Eigen::VectorXd mu = ndarray::viewAsEigen(mean);
    Eigen::MatrixXd sigma = fraction * cholesky.matrixL();
    ndarray::Array<double,1,1> parameters = ndarray::allocate(mean.getSize<0>());
    for (int i = 0; i < nComponents; ++i) {
        drawGaussian(engine, mu, sigma, parameters, workspace);
        components.push_back(
            MixtureComponent(
                1.0 / nComponents, ndarray::viewAsEigen(parameters), fraction * cholesky.matrixL()
            )
        );
    }
    return MixtureDistribution(components, dof);
}

}}}} // namespace lsst::meas::multifit::sampling
