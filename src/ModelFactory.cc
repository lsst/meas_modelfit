#include "lsst/meas/multifit/ModelFactory.h"
#include "lsst/meas/multifit/Model.h"
#include "lsst/pex/exceptions/Runtime.h"

namespace multifit = lsst::meas::multifit;

multifit::ModelFactory::RegistryMap & multifit::ModelFactory::getRegistry() {
    static RegistryMap registry;
    return registry;
}

multifit::ModelFactory::ConstPtr multifit::ModelFactory::lookup(
    std::string const & name
) {
    RegistryMap::const_iterator i = getRegistry().find(name);
    if (i == getRegistry().end()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            (boost::format("No ModelFactory associated with name '%s'.") % name).str()
        );
    }
    return i->second;
}

bool multifit::ModelFactory::declare(
    std::string const & name, 
    ModelFactory::ConstPtr const & factory
) {
    std::pair<RegistryMap::iterator,bool> result = 
        getRegistry().insert(std::make_pair(name,factory));
    return result.second;
}

multifit::Model::Ptr multifit::makeModel(
    std::string const & type, 
    ParameterVector const & linearParameters,
    ParameterVector const & nonlinearParameters
) {
    ModelFactory::ConstPtr factory = ModelFactory::lookup(type);

    return factory->makeModel(linearParameters, nonlinearParameters);
}

