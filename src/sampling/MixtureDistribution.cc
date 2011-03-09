#include "lsst/meas/multifit/sampling/MixtureDistribution.h"
#include "lsst/ndarray/eigen.h"
#include "lsst/pex/exceptions.h"

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
    double const k = 0.5 * workspace.size() * std::log(2.0 * M_PI);
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
            workspace = ndarray::viewAsEigen(record.parameters) - j->getMu();
            j->getSigma().solveTriangularInPlace(workspace);
            record.proposal += j->getNormalization() 
                * std::exp(-(k + j->_sigma.diagonal().cwise().log().sum() + 0.5 * workspace.squaredNorm()));
        }
    }
}

}}}} // namespace lsst::meas::multifit::sampling
