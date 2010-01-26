// -*- lsst-c++ -*-
/**
 * @file
 * Declaration of class ModelFactory
 */
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
 * A factory and dynamic type object for Model.
 *
 *  A ModelFactory instance provides the only means of instantiating Model
 *  objects, and also serves as a unique identifier for a particular type 
 *  of model.
 *
 *  The base ModelFactory contains a registry mapping string identifiers to
 *  ModelFactory instances.
 *
 *  As part of its role as a dynamic type object, a ModelFactory instance fixes
 *  the number of nonlinear parameters produced Models will have, and specifies
 *  a range for the number of linear parameters.
 *
 *  @sa Model
 *  @sa ModelProjection
 */
class ModelFactory : private boost::noncopyable {
public:
    typedef boost::shared_ptr<ModelFactory> Ptr;
    typedef boost::shared_ptr<ModelFactory const> ConstPtr;

    static ModelFactory::ConstPtr lookup(std::string const & name);
    
    static bool declare(
        std::string const & name, 
        ModelFactory::ConstPtr const & factory
    );

    /**
     *  Create a Model.
     *
     *  The type of model constructed, is determined by the type of factory
     *  @param linearParameters specifies the number and value of linear
     *      parameters of the model. The number of linear parameters must be
     *      within the range 
     *      [getMinLinearParameterSize():getMaxLinearParameterSize()]
     *  @param nonlinearParameters specifies the value of the nonlinear
     *      parameters of the model. The length of nonlinearParameters must
     *      equal getNonlinearParameterSize()     
     */
    virtual boost::shared_ptr<Model> makeModel(        
        ParameterVector const & linearParameters,
        ParameterVector const & nonlinearParameters
    ) const = 0;
    
    /// Number of nonlinear parameters in Models produced by this factory.
    virtual int const getNonlinearParameterSize() const = 0;

    /// Minimum number of linear parameters in Models produced by this factory.
    virtual int const getMinLinearParameterSize() const = 0;

    /// Maximum number of linear parameters in Models produced by this factory.
    virtual int const getMaxLinearParameterSize() const = 0;

    virtual ~ModelFactory() {}

private:
    typedef std::map<std::string,ModelFactory::ConstPtr> RegistryMap;
    static RegistryMap & getRegistry();
};

boost::shared_ptr<Model> makeModel(
    std::string const & type, 
    ParameterVector const & linearParameters, 
    ParameterVector const & nonlinearParameters
);


}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_MODEL_FACTORY_H
