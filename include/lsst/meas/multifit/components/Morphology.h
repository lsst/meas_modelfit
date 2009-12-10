#ifndef LSST_MEAS_MULTIFIT_COMPONENTS_MORPHOLOGY_H
#define LSST_MEAS_MULTIFIT_COMPONENTS_MORPHOLOGY_H

#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/noncopyable.hpp>

#include "lsst/afw/math/ellipses.h"
#include "lsst/afw/math/AffineTransform.h"

#include "multifit/components/MorphologyProjection.hpp"

namespace lsst {
namespace meas {
namespace multifit {

class ComponentModel;
class ComponentModelProjection;

namespace components {

/**
 *  \brief An abstract class that describes the morphology (i.e. everything but the position)
 *  of a ComponentModel.
 *
 *  \todo Better exception classes for invalid operations.
 */
class Morphology : 
    public boost::enable_shared_from_this<Morphology>, 
    private boost::noncopyable 
{
public:
    typedef boost::shared_ptr<Morphology> Ptr;
    typedef boost::shared_ptr<Morphology const> ConstPtr;

    virtual ~Morphology() {}
    
    // --- Template-mode functionality ----------------------------------------------------------------------

    /// \brief Return the minimum number of linear parameters.
    virtual int const getMinLinearParameterSize() const = 0;

    /// \brief Return the maximum number of linear parameters.
    virtual int const getMaxLinearParameterSize() const = 0;

    /// \brief Return the number of (nonlinear) morphology parameters.
    virtual int const getMorphologyParameterSize() const = 0;

    // --- In-model functionality ---------------------------------------------------------------------------

    /// \brief Return the number of linear parameters.
    int const getLinearParameterSize() const {
        if (!_linearParameterVector){
            throw LSST_EXCEPT(
                lsst::pex::exceptions::LogicErrorException,
                "Linear parameter size is indeterminate for template morphologies."
            );
        }
        return _linearParameterVector->size();
    }

    /// \brief Return a vector of the linear parameters.
    ParameterVector const & getLinearParameterVecotr() const {
        if (!_linearParameterVector){
            throw LSST_EXCEPT(
                lsst::pex::exceptions::LogicErrorException,
                "Linear parameters are not set template morphologies."
            );
        }
        return *_linearParameterVector;
    }

    /// \brief Return an ellipse core that bounds the morphology.
    virtual lsts::afw::math::ellipses::Core::Ptr computeBoundingEllipseCore() const = 0;

    /**
     *  \brief Create a new MorphologyProjection object.
     *
     *  Typically used only by the owning ComponentModel.
     */
    virtual MorphologyProjection::Ptr makeProjection(
        int kernelSize,
        lsst::afw::math::AffineTransform::ConstPtr const & transform
    ) const = 0;

protected:

    friend class multifit::ComponentModel;
    friend class multifit::ComponentModelProjection;

    /**
     *  \brief Construct a new Morphology using this as a template.
     *
     *  Typically used only by ComponentModel.
     *
     *  \param linearParameterVector The owning ComponentModel's linear parameter vector.
     *  \param morphologyParameterIter Iterator to the first Morphology-specific parameter in 
     *  the owning ComponentModel's nonlinear parameter vector. 
     */
    virtual Morphology::Ptr create(
        boost::shared_ptr<ParameterVector const> const & linearParameterVector
        ParameterConstIterator morphologyParameterIter
    ) const = 0;

    /**
     *  \brief Default-construct a Morphology object to be used as a template.
     *
     *  A public constructor that delegates to this one should be available for leaf subclasses.
     */
    Morphology() : _linearParameterVector(), _morphologyParameterIter(NULL) {}

    /**
     *  \brief Construct a Morphology object for use inside a ComponentModel.
     *
     *  \sa Morphology::create()
     */
    Morphology(
        boost::shared_ptr<ParameterVector const> const & linearParameterVector,
        ParameterConstIterator morphologyParameterIter
    ) : _linearParameterVector(linearParameterVector),
        _morphologyParameterIter(morphologyParameterIter)
    {}

    /**
     *  \brief Handle a change in the linear parameters, as propogated by the owning ComponentModel.
     */
    virtual void _handleLinearParameterChange() {}

    /**
     *  \brief Handle a change in the (nonlinear) morphology parameters, as propogated by the
     *  owning ComponentModel.
     */
    virtual void _handleMorphologyParameterChange() {}

    /// \brief Return an iterator to the Model's (nonlinear) morphology parameters.
    ParameterConstIterator _getMorphologyParameterIter() const { return _morphologyParameterIter; }

private:
    boost::shared_ptr<ParameterVector const> _linearParameterVector;
    ParameterConstIterator _morphologyParameterIter;
};

}}} // namespace lsst::meas::multifit::components

#endif // !LSST_MEAS_MULTIFIT_COMPONENTS_MORPHOLOGY_H

