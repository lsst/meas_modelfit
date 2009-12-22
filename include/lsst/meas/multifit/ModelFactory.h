#ifndef LSST_MEAS_MULTIFIT_MODEL_FACTORY_H
#define LSST_MEAS_MULTIFIT_MODEL_FACTORY_H

#include <map>
#include <string>

#include <boost/enable_shared_from_this.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/noncopyable.hpp>

#include <ndarray_fwd.hpp>

#include "lsst/meas/multifit/core.h"

namespace lsst {
namespace meas {
namespace multifit {

class Model;

/**
 *  \brief A factory and dynamic type object for Model (ABC).
 *
 *  A ModelFactory instance provides the only means of constructing Model objects,
 *  and also serves as a unique identifier for a particular type of model.
*
 *  The base ModelFactory contains a registry mapping string identifiers to
 *  ModelFactory instances.
 *
 *  As part of its role as a dynamic type object, a ModelFactory instance fixes
 *  the number of nonlinear parameters produced Models will have, and specifies
 *  a range for the number of linear parameters.
 *
 *  \sa Model
 *  \sa ModelProjection
 */
class ModelFactory : private boost::noncopyable {
public:
    typedef boost::shared_ptr<ModelFactory> Ptr;
    typedef boost::shared_ptr<ModelFactory const> ConstPtr;

    /// \brief Retrieve the factory registered with the given name.
    static ModelFactory::ConstPtr lookupFactory(std::string const & name);
    
    /**
     *  \brief Register a factory with the given name.
     *
     *  \return true on success, false if the given name is already in use.
     */
    static bool registerFactory(std::string const & name, ModelFactory::ConstPtr const & factory);

    /**
     *  \brief Create a Model.
     *
     *  \todo Add another constructor that takes direct measurement quantities for use
     *  in initializing from detection-stage results.
     */
    virtual boost::shared_ptr<Model> makeModel(
        int linearParameterSize,
        ParameterConstIterator linearParameters,
        ParameterConstIterator nonlinearParameters
    ) const = 0;
    
    /// \brief Return the number of nonlinear parameters in Models produced by this factory.
    virtual int const getNonlinearParameterSize() const = 0;

    /// \brief Return the minimum number of linear parameters in Models produced by this factory.
    virtual int const getMinLinearParameterSize() const = 0;

    /// \brief Return the maximum number of linear parameters in Models produced by this factory.
    virtual int const getMaxLinearParameterSize() const = 0;

    virtual ~ModelFactory() {}

private:
    typedef std::map<std::string,ModelFactory::ConstPtr> RegistryMap;
    static RegistryMap _registry;
};



}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_MODEL_FACTORY_H
