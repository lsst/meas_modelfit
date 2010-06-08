// -*- lsst-c++ -*-
/**
 * @file
 * Declaration of base class Morphology
 */
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
class ComponentModelProjection;

namespace components {

class FixedNonlinearMorphology;
/**
 *  An abstract class that describes the morphology (i.e. everything but the 
 *  position) of a ComponentModel.
 *
 *  Morphologies are used generally used by two classes, ComponentModel, and
 *  ComponentModelFactory. The former uses Morphology to interpret its
 *  parameters, and we refer to this mode of operation In-model mode. The latter
 *  uses Morphology as a factory for other Morphology instances, and we refer to
 *  this latter mode as Template mode
 */
class Morphology : 
    public boost::enable_shared_from_this<Morphology>, 
    private boost::noncopyable 
{
public:
    typedef boost::shared_ptr<Morphology> Ptr;
    typedef boost::shared_ptr<Morphology const> ConstPtr;

    virtual ~Morphology() {}
    
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

         
    /// Return the number of linear parameters.
    virtual int const getLinearParameterSize() const {
        return _linearParameters->size();
    }
    /// Return the number of nonlinear morphology parameters.
    virtual int const getNonlinearParameterSize() const = 0;

     
    /// Return a vector of the linear parameters.
    ParameterVector const getLinearParameters() const {
        return *_linearParameters;
    }
    ParameterVector const getNonlinearParameters() const{
        if(getNonlinearParameterSize() == 0)
            return ParameterVector();        
        return _nonlinearParameters->segment(_start, getNonlinearParameterSize());
    }

    /**
     * Return the index of the fist Nonlinear Parameters specific to this 
     * morphology
     */
    size_t const & getNonlinearParameterOffset() const {return _start;}

    /**
     * Return an iterator to the start of the morphology's nonlinear 
     * parameters.
     */
    ParameterConstIterator beginNonlinear() const { 
        return _nonlinearParameters->data() + _start; 
    }

    /**
     * Return an iterator to the end of the morphology's nonlinear parameters.
     */
    ParameterConstIterator endNonlinear() const {
        return beginNonlinear()+getNonlinearParameterSize();
    }
    
    /**
     *  Construct a new Morphology using this as a template.
     *
     *  Typically used only by ComponentModel. The parameters of this Morphology
     *  will not be shared with the new Morphology.
     *
     *  @param linearParameters The owning ComponentModel's linear 
     *      parameter vector
     *  @param nonlinearParameters The owning ComponentModel's nonlinear
     *      parameter vector
     *  @param start index of the first nonlinear parameter relevant to this
     *      morphology
     */
    virtual Morphology::Ptr create(
        boost::shared_ptr<ParameterVector const> const & linearParameters,
        boost::shared_ptr<ParameterVector const> const & nonlinearParameters,
        size_t const & start=0 
    ) const = 0;
protected:
    friend class multifit::ComponentModel;
    friend class multifit::ComponentModelProjection;
    friend class MorphologyProjection;
    friend class FixedNonlinearMorphology;
    /**
     *  Construct a Morphology object for use inside a ComponentModel.
     *
     *  @sa Morphology::create()
     */
    Morphology(
        boost::shared_ptr<ParameterVector const> const & linearParameters,
        boost::shared_ptr<ParameterVector const> const & nonlinearParameters,
        size_t const & start=0
    ) : _linearParameters(linearParameters),
        _nonlinearParameters(nonlinearParameters),
        _start(start)
    {}

    /**
     *  Handle a change in the linear parameters, as propogated by the owning 
     *  ComponentModel.
     */
    virtual void _handleLinearParameterChange() {}

    /**
     *  Handle a change in the nonlinear (morphology) parameters, as propogated
     *  by the owning ComponentModel.
     */
    virtual void _handleNonlinearParameterChange() {}

    boost::shared_ptr<ParameterVector const> _linearParameters;
    boost::shared_ptr<ParameterVector const> _nonlinearParameters;
    size_t _start;
};

}}}} // namespace lsst::meas::multifit::components

#endif // !LSST_MEAS_MULTIFIT_COMPONENTS_MORPHOLOGY_H

