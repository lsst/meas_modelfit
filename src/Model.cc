#include <ndarray.hpp>
#include <lsst/meas/multifit/Model.h>

namespace multifit = lsst::meas::multifit;

void multifit::Model::getLinearParameters(
        ParameterIterator parameters
) const {
    Eigen::VectorXd const & vector = getLinearParameters();
    std::copy(vector.data(),vector.data()+vector.size(),parameters);
}

void multifit::Model::getNonlinearParameters(
        ParameterIterator parameters
) const {
    Eigen::VectorXd const & vector = getNonlinearParameters();
    std::copy(vector.data(),vector.data()+vector.size(),parameters);
}

multifit::Model::Definition::Definition(
        Factory::ConstPtr const & factory,
        ParameterConstIterator linearParameterBegin,
        ParameterConstIterator const linearParameterEnd,
        ParameterConstIterator nonlinearParameterBegin
) : _factory(factory),
    _linearParameters(linearParameterEnd-linearParameterBegin),
    _nonlinearParameters(factory->getNonlinearParameterSize())
{
    setLinearParameters(linearParameterBegin);
    setNonlinearParameters(nonlinearParameterBegin);
}

multifit::Model::Model(
        Definition const & definition,
        int activeProducts,
        WcsConstPtr const &wcs,
        FootprintConstPtr const &footprint,
        double photFactor
) : _activeProducts(activeProducts),
    _photFactor(photFactor),
    _definition(new Definition(definition)),
    _footprint(footprint),
    _wcs(wcs)
{}

multifit::Model::Model(Model const & other) : 
    _activeProducts(other._activeProducts),
    _photFactor(other._photFactor),
    _definition(new Definition(*other._definition)),
    _footprint(other._footprint),
    _wcs(other._wcs)
{}
