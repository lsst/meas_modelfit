#ifndef LSST_MEAS_MULTIFIT_COMPONENTS_MORPHOLOGY_H
#define LSST_MEAS_MULTIFIT_COMPONENTS_MORPHOLOGY_H

#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/noncopyable.hpp>

#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/components/MorphologyProjection.h"

namespace lsst {
namespace meas {
namespace multifit {

class ComponentModel;
class ComponentModelFactory;
class ComponentModelProjection;

namespace components {

/**
 *  An abstract class that describes the morphology (i.e. everything but the position)
 *  of a ComponentModel.
 *
 *  Morphologies are used generally used by two classes, ComponentModel, and
 *  ComponentModelFactory. The former uses Morphology to interpret its
 *  parameters, and we refer to this mode of operation In-model mode. The latter uses 
 *  Morphology as a factory for other Morphology instances, and we refer to this latter
 *  mode as Template mode
 */
class Morphology : 
    public boost::enable_shared_from_this<Morphology>, 
    private boost::noncopyable 
{
public:
    typedef boost::shared_ptr<Morphology> Ptr;
    typedef boost::shared_ptr<Morphology const> ConstPtr;

    virtual ~Morphology() {}
    


    /// Return the number of linear parameters.
    int const getLinearParameterSize() const {
        if (!_linearParameterVector){
            throw LSST_EXCEPT(
                lsst::pex::exceptions::LogicErrorException,
                "Linear parameter size is indeterminate for template morphologies."
            );
        }
        return _linearParameterVector->size();
    }

    /// Return a vector of the linear parameters.
    ParameterVector const & getLinearParameterVector() const {
        if (!_linearParameterVector){
            throw LSST_EXCEPT(
                lsst::pex::exceptions::LogicErrorException,
                "Linear parameters are not set template morphologies."
            );
        }
        return *_linearParameterVector;
    }
    /**
     * @name In-model Mode Functionality
     *
     * These methods are useful to ComponentModel, which delegates to a
     * Morpholgy. 
     */
    //@{
    
    /**
     * Return an ellipse core that bounds the morphology.
     *
     * @sa ComponentModel::computeBoundingEllipse
     */
    virtual lsst::afw::geom::ellipses::Core::Ptr computeBoundingEllipseCore() const = 0;

    /**
     *  Create a new MorphologyProjection object.
     *
     *  @sa ComponentModel::makeProjection
     */
    virtual MorphologyProjection::Ptr makeProjection(
        lsst::afw::geom::Extent2I const & kernelSize,
        lsst::afw::geom::AffineTransform::ConstPtr const & transform
    ) const = 0;
    //@}
protected:

    friend class multifit::ComponentModel;
    friend class multifit::ComponentModelFactory;
    friend class multifit::ComponentModelProjection;

    /**
     *  Construct a new Morphology using this as a template.
     *
     *  Typically used only by ComponentModel. The parameters of this Morphology
     *  will not be shared with the new Morphology.
     *
     *  @param linearParameterVector The owning ComponentModel's linear parameter vector
     *  @param nonlinearParameterIter Iterator to the first Morphology-specific
     *   parameter in the owning ComponentModel's nonlinear parameter vector.
     */
    virtual Morphology::Ptr create(
        boost::shared_ptr<ParameterVector const> const & linearParameterVector,
        ParameterConstIterator nonlinearParameterIter
    ) const = 0;

    /**
     *  Default-construct a Morphology object to be used as a template.
     *
     *  A public constructor that delegates to this one should be available for leaf subclasses.
     */
    Morphology() : _linearParameterVector(), _nonlinearParameterIter(NULL) {}

    /**
     *  Construct a Morphology object for use inside a ComponentModel.
     *
     *  @sa Morphology::create()
     */
    Morphology(
        boost::shared_ptr<ParameterVector const> const & linearParameterVector,
        ParameterConstIterator nonlinearParameterIter
    ) : _linearParameterVector(linearParameterVector),
        _nonlinearParameterIter(nonlinearParameterIter)
    {}

    /**
     * @name Template Mode Functionality
     *
     * These methods are only useful when a morphology is being used as a
     * template to construct other Morphology instances. Generally, this occurs
     * within a ComponentModelFactory.
     * 
     * Morphologies may be able to interpret a range of linear paramter sizes.
     * By allowing the morpholgy, rather than the model, to dictate this value,
     * the components of a ComponentModel encapsulate all the bhevaiour of the
     * model. This allows us to derive new model types simply by deriving new
     * morphologies, rather than whole new models.
     */
    //@{
    /// Return the minimum number of linear parameters.
    virtual int const getMinLinearParameterSize() const = 0;

    /// Return the maximum number of linear parameters.
    virtual int const getMaxLinearParameterSize() const = 0;

    /// Return the number of (nonlinear) morphology parameters.
    virtual int const getNonlinearParameterSize() const = 0;
    //@}

    /**
     *  Handle a change in the linear parameters, as propogated by the owning ComponentModel.
     */
    virtual void _handleLinearParameterChange() {}

    /**
     *  Handle a change in the nonlinear (morphology) parameters, as propogated by the
     *  owning ComponentModel.
     */
    virtual void _handleNonlinearParameterChange() {}

    /// Return an iterator to the Model's (nonlinear) morphology parameters.
    ParameterConstIterator _getNonlinearParameterIter() const { return _nonlinearParameterIter; }

private:
    boost::shared_ptr<ParameterVector const> _linearParameterVector;
    ParameterConstIterator _nonlinearParameterIter;

};

}}}} // namespace lsst::meas::multifit::components

#endif // !LSST_MEAS_MULTIFIT_COMPONENTS_MORPHOLOGY_H

