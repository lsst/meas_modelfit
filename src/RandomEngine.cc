#include "lsst/meas/multifit/RandomEngine.h"

namespace lsst { namespace meas { namespace multifit {

RandomEngine::RandomEngine() :
    _engine(),
    _normal(_engine, boost::normal_distribution<double>(0.0, 1.0)),
    _uniform(_engine, boost::uniform_real<double>(0.0, 1.0))
{}

RandomEngine::RandomEngine(RandomEngine const & other) :
    _engine(other._engine),
    _normal(_engine, boost::normal_distribution<double>(0.0, 1.0)),
    _uniform(_engine, boost::uniform_real<double>(0.0, 1.0))
{}

RandomEngine & RandomEngine::operator=(RandomEngine const & other) {
    if (&other != this) {
        _engine = other._engine;
    }
    return *this;
}

}}} // namespace lsst::meas::multifit
