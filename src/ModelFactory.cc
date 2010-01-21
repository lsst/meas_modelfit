#include "lsst/meas/multifit/ModelFactory.h"
#include "lsst/meas/multifit/Model.h"
#include "lsst/pex/exceptions/Runtime.h"

namespace multifit = lsst::meas::multifit;

/**
 * Get a dictionary of ModelFactories indexed by name
 */
multifit::ModelFactory::RegistryMap & multifit::ModelFactory::getRegistry() {
    static RegistryMap registry;
    return registry;
}

/**
 * Retrieve the ModelFactory of the given name
 *
 * @param name desired ModelFactory.
 * @throw lsst::pex::exceptions::NotFoundException if no factory of that name
 * exists
 */
multifit::ModelFactory::ConstPtr multifit::ModelFactory::lookup(
    std::string const & name
) {
    RegistryMap::const_iterator i = getRegistry().find(name);
    if (i == getRegistry().end()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::NotFoundException,
            (boost::format("No ModelFactory associated with name '%s'.") % name).str()
        );
    }
    return i->second;
}

/**
 * Add a ModelFactory to the registry
 *
 * When choosing a name for a ModelFactory, used the type of Model they will be
 * used to construct. 
 * e.g. "PointSource" names a factory used to construct Point Source Model.
 *
 * @param name factory will be added to the registry of \c ModelFactory objects
 * under this name
 * @param factory the \c ModelFactory to add to the registry.
 * @return true if declaration was successful 
 */
bool multifit::ModelFactory::declare(
    std::string const & name, 
    ModelFactory::ConstPtr const & factory
) {
    std::pair<RegistryMap::iterator,bool> result = 
        getRegistry().insert(std::make_pair(name,factory));
    return result.second;
}

/**
 * Construct a Model using the ModelFactory of declared by name type
 *
 * @throw lsst::pex::exceptions::NotFoundException if no ModelFactory has been
 * declared with given name
 */
multifit::Model::Ptr multifit::makeModel(
    std::string const & type, 
    ParameterVector const & linearParameters,
    ParameterVector const & nonlinearParameters
) {
    ModelFactory::ConstPtr factory = ModelFactory::lookup(type);

    return factory->makeModel(linearParameters, nonlinearParameters);
}

