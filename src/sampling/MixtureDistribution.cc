#include "lsst/meas/multifit/sampling/MixtureDistribution.h"
#include "lsst/ndarray/eigen.h"
#include "lsst/pex/exceptions.h"

#include <Eigen/Array>
#include <Eigen/Cholesky>

#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>

namespace lsst { namespace meas { namespace multifit { namespace sampling {

namespace {

double evaluateGaussian(
    MixtureComponent const & component,
    ndarray::Array<double const,1,1> const & parameters, 
    Eigen::VectorXd & workspace
) {
    static double const k = 0.5 * std::log(2.0 * M_PI);
    workspace = ndarray::viewAsEigen(parameters) - component.getMu();
    component.getSigma().part<Eigen::LowerTriangular>().solveTriangularInPlace(workspace);
    return component.getNormalization() * std::exp(
        -k
        - workspace.size() 
        - 0.5 * workspace.squaredNorm() 
        - component.getSigma().diagonal().cwise().log().sum()
    );
}

} // anonymous

MixtureDistribution::MixtureDistribution(ComponentList const & components) : _components(components) {
    double total = 0.0;
    for (ComponentList::iterator j = _components.begin(); j != _components.end(); ++j) {
        total += j->getNormalization();
    }
    for (ComponentList::iterator j = _components.begin(); j != _components.end(); ++j) {
        j->_normalization /= total;
    }
}

void MixtureDistribution::draw(Table const & table, RandomEngine & engine) const {
    if (table.getParameterSize() != getParameterSize()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LengthErrorException,
            "Number of parameters in table does not match number of parameters in distribution"
        );
    }
    typedef boost::variate_generator< RandomEngine &, boost::uniform_real<double> > UniformGenerator;
    typedef boost::variate_generator< RandomEngine &, boost::normal_distribution<double> > NormalGenerator;
    UniformGenerator uniform(engine, boost::uniform_real<double>(0.0, 1.0));
    NormalGenerator normal(engine, boost::normal_distribution<double>(0.0, 1.0));
    Eigen::VectorXd workspace(getParameterSize());
    table.clear();
    for (int i = 0; i < table.getTableSize(); ++i) {
        Record record = table[i];

        // First draw the parameter point
        for (int n = 0; n < workspace.size(); ++n) {
            workspace[n] = normal();
        }
        double u = uniform();
        double v = 0.0;
        ComponentList::const_iterator j = _components.begin();
        for (; j != _components.end(); ++j) {
            v += j->getNormalization();
            if (v >= u) break;
        }
        if (j == _components.end()) --j; // only needed because of round-off error
        ndarray::viewAsEigen(record.parameters) = j->getMu() + j->getSigma() * workspace;

        // Now evaluate the probability at that point.
        for (j = _components.begin(); j != _components.end(); ++j) {
            record.importance += evaluateGaussian(*j, record.parameters, workspace);
        }
    }
}

void MixtureDistribution::update(ConstTable const & table) {
    Eigen::VectorXd rho(table.getTableSize());
    Eigen::VectorXd tau(table.getTableSize());
    Eigen::VectorXd mu(table.getParameterSize());
    Eigen::MatrixXd sigma(table.getParameterSize(), table.getParameterSize());
    Eigen::VectorXd workspace(table.getParameterSize());
    for (ComponentList::iterator j = _components.begin(); j != _components.end(); ++j) {
        double alpha = 0.0;
        for (int i = 0; i < table.getTableSize(); ++i) {
            ConstRecord record = table[i];
            rho[i] = evaluateGaussian(*j, record.parameters, workspace) / record.importance;
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
        j->_normalization = alpha;
        j->_mu = mu;
        Eigen::LLT<Eigen::MatrixXd> cholesky(sigma);
        j->_sigma = cholesky.matrixL();
    }
}

}}}} // namespace lsst::meas::multifit::sampling
