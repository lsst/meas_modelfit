#include "lsst/meas/multifit/sampling/MixtureDistribution.h"
#include "lsst/ndarray/eigen.h"

#include <Eigen/Array>

#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>

namespace lsst { namespace meas { namespace multifit { namespace sampling {

MixtureDistribution::MixtureDistribution(ComponentList const & components) : _components(components) {
    double total = 0.0;
    for (ComponentList::iterator j = _components.begin(); j != _components.end(); ++j) {
        total += j->getNormalization();
    }
    for (ComponentList::iterator j = _components.begin(); j != _components.end(); ++j) {
        j->_normalization /= total;
    }
}

void MixtureDistribution::draw(
    lsst::ndarray::Array<double,2,2> const & points,
    RandomEngine & engine
) const {
    typedef boost::variate_generator< RandomEngine &, boost::uniform_real<double> > UniformGenerator;
    typedef boost::variate_generator< RandomEngine &, boost::normal_distribution<double> > NormalGenerator;
    UniformGenerator uniform(engine, boost::uniform_real<double>(0.0, 1.0));
    NormalGenerator normal(engine, boost::normal_distribution<double>(0.0, 1.0));
    Eigen::VectorXd workspace(points.getSize<1>());
    for (lsst::ndarray::Array<double,2,2>::Iterator i = points.begin(); i != points.end(); ++i) {
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
        ndarray::viewAsEigen(*i) = j->getMu() + j->getSigma() * workspace;
    }
}

void MixtureDistribution::evaluate(
    lsst::ndarray::Array<double,1,1> const & probability,
    lsst::ndarray::Array<double const,2,2> const & points
) const {
    Eigen::MatrixXd workspace(points.getSize<1>(), points.getSize<0>());
    ndarray::EigenView<double,1,1> p(probability);
    p.setZero();
    Eigen::RowVectorXd::ConstantReturnType ones = Eigen::RowVectorXd::Ones(points.getSize<0>());
    double const k1 = 0.5 * points.getSize<1>() * std::log(2.0 * M_PI);
    for (ComponentList::const_iterator j = _components.begin(); j != _components.end(); ++j) {
        workspace = ndarray::viewAsTransposedEigen(points) - j->getMu() * ones;
        j->getSigma().solveTriangularInPlace(workspace);
        double k2 = k1 + j->_sigma.diagonal().cwise().log().sum();
        for (int n = 0; n < p.size(); ++n) {
            p[n] += j->getNormalization() * std::exp(-(k2 + 0.5 * workspace.col(n).squaredNorm()));
        }
    }
}


}}}} // namespace lsst::meas::multifit::sampling
