#ifndef LSST_MEAS_MULTIFIT_MODEL_PROJECTION_H
#define LSST_MEAS_MULTIFIT_MODEL_PROJECTION_H

#include <Eigen/Core>
#include <ndarray.hpp>

#include "lsst/meas/multifit/core.h"
#include "lsst/meas/multifit/Model.h"

namespace lsst {
namespace meas {
namespace multifit {

/**
 *  \brief A projection of a Model to a particular set of observing conditions 
 *  (ABC).
 *
 *  ModelProjection represents a Model evaluated on a particular exposure, and
 *  provides functionality to create an image representation of the Model.
 *  Various derivatives of the image representation with respect to different
 *  model parameters and calibration parameters can also be computed, with the
 *  necessary change-of-variables terms automatically included.
 *  
 *  \sa Model
 *  \sa ModelFactory
 */
class ModelProjection : private boost::noncopyable {
public:

    typedef boost::shared_ptr<ModelProjection> Ptr;
    typedef boost::shared_ptr<ModelProjection const> ConstPtr;

    virtual ~ModelProjection() {}

    /// \brief Return the Model instance this is a projection of.
    Model::ConstPtr getModel() const { return _model; }

    /// \brief Return the number of linear parameters.
    int const getLinearParameterSize() const { return _model->getLinearParameterSize(); }

    /// \brief Return the number of nonlinear parameters.
    int const getNonlinearParameterSize() const { return _model->getNonlinearParameterSize(); }

    /// \brief Return the number of parameters that specify the coordinate transformation.
    virtual int const getWcsParameterSize() const = 0;

    /// \brief Return the number of parameters that specify the PSF.
    virtual int const getPsfParameterSize() const = 0;

    /**
     *  @name PublicProductComputers
     *
     *  Each of these computes the Footprint-compressed product as a row-major array
     *  with the last dimension corresponding to the Footprint pixel index and the first
     *  dimension (if any) corresponding to the parameter index.
     */
    //@{
    ndarray::Array<Pixel const,1,1> computeModelImage();
    ndarray::Array<Pixel const,2,1> computeLinearParameterDerivative();
    ndarray::Array<Pixel const,2,1> computeNonlinearParameterDerivative();
    ndarray::Array<Pixel const,2,1> computeWcsParameterDerivative();
    ndarray::Array<Pixel const,2,1> computePsfParameterDerivative();
    //@}

    /**
     *  @name ProductBufferSetters
     *
     *  The compute[Product]() member functions can be set to fill and return externally
     *  allocated (and possibly shared) data.  After setting a buffer, only the
     *  ModelProjection should be allowed to modify it.
     */
    //@{
    void setModelImageBuffer(ndarray::Array<Pixel,1,1> const & buffer);
    void setLinearParameterDerivativeBuffer(ndarray::Array<Pixel,2,1> const & buffer);
    void setNonlinearParameterDerivativeBuffer(ndarray::Array<Pixel,2,1> const & buffer);
    void setWcsParameterDerivativeBuffer(ndarray::Array<Pixel,2,1> const & buffer);
    void setPsfParameterDerivativeBuffer(ndarray::Array<Pixel,2,1> const & buffer);
    //@}

    /// \brief Return the World Coordinate System object for this projection.
    WcsConstPtr const & getWcs() const { return _wcs; }

    /// \brief Return the Footprint the image representation will be computed on.
    FootprintConstPtr const & getFootprint() const { return _footprint; }

protected:


    /// \brief Construct a projection.
    ModelProjection(
        Model::ConstPtr const & model,
        WcsConstPtr const & wcs,
        FootprintConstPtr const & footprint
    );

    /**
     *  @name ProtectedProductComputers
     *
     *  Each of these computes a Footprint-compressed derivative of the projected Model's
     *  image representation with respect to a different set of parameters.  The base
     *  ModelProjection guarantees it will call these with appropriately-sized arrays,
     *  and only when the associated product flag is enabled.
     */
    //@{
    virtual void _computeModelImage(ndarray::Array<Pixel,1,1> const & vector);
    virtual void _computeLinearParameterDerivative(ndarray::Array<Pixel,2,1> const & matrix) = 0;
    virtual void _computeNonlinearParameterDerivative(ndarray::Array<Pixel,2,1> const & matrix) = 0;
    virtual void _computeWcsParameterDerivative(ndarray::Array<Pixel,2,1> const & matrix) = 0;
    virtual void _computePsfParameterDerivative(ndarray::Array<Pixel,2,1> const & matrix) = 0;
    //@}

    /**
     *  \brief Handle a linear parameter change broadcast from the associated Model.
     *
     *  Subclasses which override must call the base class implementation.
     */
    virtual void _handleLinearParameterChange();

    /**
     *  \brief Handle a nonlinear parameter change broadcast from the associated Model.
     *
     *  Subclasses which override must call the base class implementation.
     */
    virtual void _handleNonlinearParameterChange();

private:
    friend class Model;
    
    enum ProductFlag {
        MODEL_IMAGE = 1<<0,
        LINEAR_PARAMETER_DERIVATIVE = 1<<1,
        NONLINEAR_PARAMETER_DERIVATIVE = 1<<2,
        WCS_PARAMETER_DERIVATIVE = 1<<3,
        PSF_PARAMETER_DERIVATIVE = 1<<4,
    };

    int _validProducts;

    Model::ConstPtr _model;
    FootprintConstPtr _footprint;
    WcsConstPtr _wcs;

    ndarray::Array<Pixel, 1, 1> _modelImage;
    ndarray::Array<Pixel, 2, 1> _linearParameterDerivative;
    ndarray::Array<Pixel, 2, 1> _nonlinearParameterDerivative;
    ndarray::Array<Pixel, 2, 1> _wcsParameterDerivative;
    ndarray::Array<Pixel, 2, 1> _psfParameterDerivative;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_MODEL_PROJECTION_H
