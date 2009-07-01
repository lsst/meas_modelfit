#include "lsst/meas/multifit/Model.h"

///////////////////////////////////////////////////////////////////////////////
// Factory methods

namespace multifit = lsst::meas::multifit;

multifit::ModelFactoryBase::ModelFactoryBase(
        std::string const & type, 
        std::string const & constraint) {
    RegistryMap& registry = getRegistry();
    registry.insert(std::make_pair(std::make_pair(type,constraint),this));
}

multifit::ModelFactoryBase::RegistryMap& multifit::ModelFactoryBase::getRegistry() {
    static RegistryMap instance;
    return instance;
}

multifit::ObjectModel* multifit::createObjectModel(
        std::string const & type, 
        std::string const & constraint,
        ndarray::ArrayRef<ObjectModel::Real,1,1> const & nonlinear_params, 
        ndarray::ArrayRef<ObjectModel::Real,1,1> const & linear_params) {
    ModelFactoryBase::RegistryMap& registry = ModelFactoryBase::getRegistry();
    ModelFactoryBase::RegistryIterator iter = 
            registry.find(std::make_pair(type,constraint));
    if(iter == registry.end()) 
    {
        throw std::runtime_error("Model of type \"" + type + 
                "\" with constraint \"" + constraint + "\" is not implemented");
    }
    return iter->second->create(nonlinear_params,linear_params);
}

