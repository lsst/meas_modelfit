#ifndef LSST_MEAS_MULTIFIT_MODEL_H
#define LSST_MEAS_MULTIFIT_MODEL_H

#include <Eigen/Core>

#include <list>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/make_shared.hpp>
#include <boost/weak_ptr.hpp>

#include <ndarray_fwd.hpp>

#include <lsst/meas/multifit/core.h>

namespace lsst{
namespace meas {
namespace multifit {

namespace projections {

class ModelProjection;

} // namespace projections

class ModelEvaluator;
class ModelFactory;

/**
 *  \brief A model for an astronomical object (ABC).
 *
 *  A Model's parameters are always in "global" (typically celestial) coordinates.
 *  A projection of the Model to a particular image coordinate system (along with
 *  convolution by the appropriate PSF and application of other observational effects)
 *  is represented by an instance of ModelProjection, which will generally be subclassed in
 *  tandem with Model and ModelFactory.  Model and ModelProjection participate in a
 *  Observer design pattern, with a Model broadcasting changes in its parameters to its
 *  associated ModelProjection objects.
 *
 *  \sa ModelFactory
 *  \sa ModelProjection
 */
class Model : protected boost::enable_shared_from_this<Model> {
public:
    typedef boost::shared_ptr<Model> Ptr;
    typedef boost::shared_ptr<Model const> ConstPtr;

    /**
     *  \brief Create a Footprint that would contain a projection of the Model.
     */
    virtual Footprint::Ptr computeProjectionFootprint(
        Kernel::ConstPtr const & kernel,
        Wcs::ConstPtr const & wcs,
        double photFactor
    ) const = 0;

    /**
     *  \brief Create an image-coordinate bounding box that would contain a
     *  projection of the Model.
     */
    virtual lsst::afw::image::BBox computeProjectionEnvelope(
        Kernel::ConstPtr const & kernel,
        Wcs::ConstPtr const & wcs,
        double photFactor
    ) const = 0;

    /**
     *  \brief Create an ra/dec bounding ellipse for the Model.
     */
    virtual lsst::afw::math::ellipses::Ellipse::Ptr computeBoundingEllipse() const = 0;

    /// \brief Return a vector of the linear parameters.
    ParameterVector const & getLinearParameters() const { 
        return *_linearParameterVector; 
    }

    /// \brief Return an iterator to the beginning of the linear parameters.
    ParameterConstIterator getLinearParameterIter() const { 
        return _linearParameterVector->data(); 
    }

    /// \brief Return the number of linear parameters.
    int const getLinearParameterSize() const { 
        return _linearParameterVector->size(); 
    }

    /// \brief Set the linear parameters.
    void setLinearParameters(ParameterConstIterator parameterIter);

    /// \brief Return a vector of the nonlinear parameters.
    ParameterVector const & getNonlinearParameters() const { 
        return *_nonlinearParameterVector; 
    }

    /// \brief Return an iterator to the beginning of the nonlinear parameters.
    ParameterConstIterator getNonlinearParameterIter() const { 
        return _nonlinearParameterVector->data(); 
    }

    /// \brief Return the number of nonlinear parameters.
    int const getNonlinearParameterSize() const { 
        return _nonlinearParameterVector->size(); 
    }

    /// \brief Set the nonlinear parameters.
    void setNonlinearParameters(ParameterConstIterator parameterIter);

    /**
     *  \brief Create a new Model with the same type and parameters.
     *
     *  Associated ModelProjection objects are not copied or shared;
     *  the new Model will not have any associated ModelProjections.
     */
    virtual Model::Ptr clone() const = 0;

    virtual ~Model() {}

protected:

    /** 
     *  \brief Create a ModelProjection object associated with this.
     *
     *  All public ModelProjection creation should use this function (or one that delegates to it).
     *
     *  makeProjection() delegates to _createProjection() the actual construction of a ModelProjection
     *  of the appropriate subclass, then:
     *  - Initializes the data members (model, WCS, footprint, phot. factor) of the new ModelProjection.
     *  - Calls ModelProjection::_setKernel(kernel).
     *  - Calls ModelProjection::enableProducts(activeProducts).
     *  - Calls ModelProjection::_handleLinearParameterChange().
     *  - Calls ModelProjection::_handleNonlinearParameterChange().
     *  - Registers the ModelProjection as an observer of this.
     */
    virtual boost::shared_ptr<projections::ModelProjection> makeProjection(
        Kernel::ConstPtr const & kernel,
        Wcs::ConstPtr const & wcs,
        Footprint::ConstPtr const & footprint,
        double photFactor,
        int activeProducts = 0
    ) const = 0;

    /// \brief Notify all associated ModelProjections that the linear parameters have changed.
    void _broadcastLinearParameterChange() const;

    /// \brief Notify all associated ModelProjections that the nonlinear parameters have changed.
    void _broadcastNonlinearParameterChange() const;

    /**
     *  \brief Provide additional code for setLinearParameters().
     *
     *  This will be called by setLinearParameters, after the parameter vector has been updated and
     *  before the call to _broadcastLinearParameterChange().
     */
    virtual void _handleLinearParameterChange() {}

    /**
     *  \brief Provide additional code for setNonlinearParameters().
     *
     *  This will be called by setNonlinearParameters, after the parameter vector has been updated
     *  and before the call to _broadcastNonlinearParameterChange().
     */
    virtual void _handleNonlinearParameterChange() {}

    /// \brief Initialize the Model and allocate space for the parameter vectors.
    Model(int linearParameterSize, int nonlinearParamterSize);

    /**
     * \brief Deep-copy the Model.
     *
     * This is a deep copy of the model parameters, projections will not be
     * associated with the new Model
     */ 
    explicit Model(Model const & model);

    boost::shared_ptr<ParameterVector> _linearParameterVector;
    boost::shared_ptr<ParameterVector> _nonlinearParameterVector;

private:
    typedef boost::weak_ptr<projections::ModelProjection> ProjectionWeakPtr;
    typedef std::list<ProjectionWeakPtr> ProjectionList;

    friend class ModelFactory;
    friend class ModelEvaluator;

    void operator=(Model const & other) { assert(false); } // Assignment disabled.

    boost::shared_ptr<ModelFactory const> _factory;
    mutable ProjectionList _projectionList;
};

} // namespace multifit

#endif // !LSST_MEAS_MULTIFIT_MODEL_H
