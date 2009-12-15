#ifndef LSST_MEAS_MULTIFIT_COMPONENT_MODEL_FACTORY_H
#define LSST_MEAS_MULTIFIT_COMPONENT_MODEL_FACTORY_H

#include <boost/make_shared.hpp>

#include "lsst/meas/multifit/Model.h"
#include "lsst/meas/multifit/ModelFactory.h"
#include "lsst/meas/multifit/components/Astrometry.h"
#include "lsst/meas/multifit/components/Morphology.h"

namespace lsst { 
namespace meas {
namespace multifit {

class ComponentModel;
 
/**
 *  \brief A factory and dynamic type object for ComponentModel (ABC).
 *
 *  \sa ComponentModel
 *  \sa ComponentModelProjection
 *
 *  \todo better documentation
 */
class ComponentModelFactory : public ModelFactory {
public:

    typedef boost::shared_ptr<ComponentModelFactory> Ptr;
    typedef boost::shared_ptr<ComponentModelFactory const> ConstPtr;

    /**
     *  \brief Construct a new ComponentModelFactory from template Astrometry 
     *  and Morphology objects.
     */
    static ModelFactory::ConstPtr create(
        components::Astrometry::ConstPtr const & astrometryTemplate,
        components::Morphology::ConstPtr const & morphologyTemplate
    ) {
        return ModelFactory::ConstPtr(
            new ComponentModelFactory(astrometryTemplate,morphologyTemplate)
        );
    }
    
    /**
     *  \brief Create a ComponentModel.
     */
    virtual Model::Ptr makeModel(
        int linearParameterSize,
        ParameterConstIterator linearParameterIter,
        ParameterConstIterator nonlinearParameterIter
    ) const;

    /**
     * \brief compute the number of nonlinear parameters in Models produced by 
     * this factory.
     */
    virtual int const getNonlinearParameterSize() const {
        return _astrometryTemplate->getAstrometryParameterSize() 
            + _morphologyTemplate->getMorphologyParameterSize();
    }

    /**
     * \brief Return the minimum number of linear parameters in Models produced 
     * by this factory.
     */
    virtual int const getMinLinearParameterSize() const {
        return _morphologyTemplate->getMinLinearParameterSize();
    }

    /**
     * \brief Return the maximum number of linear parameters in Models produced 
     * by this factory.
     */
     virtual int const getMaxLinearParameterSize() const {
        return _morphologyTemplate->getMaxLinearParameterSize();
    }

    /**
     * \brief Return the Astrometry used as a template in constructing 
     * ComponentModels.
     */
    components::Astrometry::ConstPtr getAstrometryTemplate() const { 
        return _astrometryTemplate; 
    }

    /**
     * \brief Return the Morphology used as a template in constructing 
     * ComponentModels.
     */
    components::Morphology::ConstPtr getMorphologyTemplate() const { 
        return _morphologyTemplate; 
    }
    
    virtual ~ComponentModelFactory() {}

private:
    ComponentModelFactory(
        components::Astrometry::ConstPtr const & astrometryTemplate,
        components::Morphology::ConstPtr const & morphologyTemplate
    ) : _astrometryTemplate(astrometryTemplate), 
        _morphologyTemplate(morphologyTemplate) 
    {}

    components::Astrometry::ConstPtr _astrometryTemplate;
    components::Morphology::ConstPtr _morphologyTemplate;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_COMPONENT_MODEL_FACTORY_H
