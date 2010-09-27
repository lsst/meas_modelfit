// -*- lsst-c++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
/**
 * @file
 * Declaration of class ModelProjection.
 */
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
 *  A projection of a Model to a particular set of observing conditions 
 *
 *  ModelProjection represents a Model evaluated on a particular exposure, and
 *  provides functionality to create an image representation of the Model.
 *  Various derivatives of the image representation with respect to different
 *  model parameters and calibration parameters can also be computed, with the
 *  necessary change-of-variables terms automatically included.
 *  
 *  @sa Model
 *  @sa ModelFactory
 */
class ModelProjection : private boost::noncopyable {
public:

    typedef boost::shared_ptr<ModelProjection> Ptr;
    typedef boost::shared_ptr<ModelProjection const> ConstPtr;

    virtual ~ModelProjection() {}

    /// Model instance this is a projection of.
    Model::ConstPtr getModel() const { return _model; }

    /// Number of linear parameters.
    int const getLinearParameterSize() const { 
        return _model->getLinearParameterSize(); 
    }

    /// Number of nonlinear parameters.
    int const getNonlinearParameterSize() const { return _model->getNonlinearParameterSize(); }

    /**
     *  @name Public Product Computers
     *
     *  Each of these functions compute a footprint-compressed product as a 
     *  row-major array with the last dimension corresponding to the 
     *  footprint-mapped pixel index and the first dimension (if any) 
     *  corresponding to the parameter index.
     *
     *  Derived classes of Model should not reimplement these methods, instead
     *  they should implement the protected product computers, which these
     *  functions delegate to.
     *
     *  @sa _computeModelImage
     *  @sa _computeLinearParameterDerivative
     *  @sa _computeNonlinearParameterDerivative
     */
    //@{
#ifndef SWIG
    ndarray::Array<Pixel const,1,1> computeModelImage();
    ndarray::Array<Pixel const,2,1> computeLinearParameterDerivative();
    ndarray::Array<Pixel const,2,1> computeNonlinearParameterDerivative();
#endif
    //@}

    /**
     *  @name Product Enabled Checkers
     *
     *  These functions allow users to check if this model supports the various
     *  products. Products will only be computed when the repective checker
     *  returns true
     */
    //@{
    bool hasModelImage() const { return true; }
    bool hasLinearParameterDerivative() const { return getLinearParameterSize() > 0; }
    bool hasNonlinearParameterDerivative() const { return getNonlinearParameterSize() > 0; }
    //@}

    /// The WCS associated with the observation this projection represents.
    lsst::afw::image::Wcs::ConstPtr const & getWcs() const { return _wcs; }


    
    /**
     * Footprint the image representation will be computed on.
     *
     * This footprint determines the size of all computed products, and the
     * mapping from compressed to uncompressed images
     */
    boost::shared_ptr<lsst::afw::detection::Footprint const> const & getFootprint() const { 
        return _footprint; 
    }

protected:
    ModelProjection(
        Model::ConstPtr const & model,
        lsst::afw::geom::AffineTransform const & transform,
        CONST_PTR(lsst::afw::detection::Footprint) const & footprint
    );

    ModelProjection(
        Model::ConstPtr const & model,
        CONST_PTR(lsst::afw::image::Wcs) const & wcs,
        CONST_PTR(lsst::afw::detection::Footprint) const & footprint
    );

    /**
     * The AffineTransfrom associated with the observation this projection 
     * represents
     */
    lsst::afw::geom::AffineTransform const & getTransform() const {
        return _skyToPixelTransform;
    }

    /**
     *  @name Protected Product Computers
     *
     *  Each of these computes a Footprint-compressed derivative of the 
     *  projected Model's image representation with respect to a different set 
     *  of parameters.  The base ModelProjection guarantees it will call these 
     *  with appropriately-sized arrays, and only when the associated product 
     *  flag is enabled.
     */
    //@{
    virtual void _computeModelImage(ndarray::Array<Pixel,1,1> const & vector);
    virtual void _computeLinearParameterDerivative(
        ndarray::Array<Pixel,2,1> const & matrix
    ) = 0;
    virtual void _computeNonlinearParameterDerivative(
        ndarray::Array<Pixel,2,1> const & matrix
    ) = 0;

    //@}

    /**
     *  Handle a linear parameter change broadcast from the associated Model.
     *
     *  Subclasses which override must call the base class implementation.
     */
    virtual void _handleLinearParameterChange();

    /**
     *  Handle a nonlinear parameter change broadcast from the associated Model.
     *
     *  Subclasses which override must call the base class implementation.
     */
    virtual void _handleNonlinearParameterChange();

private:
    friend class Model;
    friend class ModelEvaluator;

    /**
     *  @name Product Buffer Setters
     *
     *  The memory buffers on which the compute[Product]() member functions 
     *  operate can be set externally. This allows ModelEvaluator to assign 
     *  blocks of larger buffer to a ModelProjection. ModelProjection takes
     *  ownership of the buffer given to it, and the buffer should not be
     *  modified externally from the point on.
     *
     *  @sa ModelEvaluator
     */
    //@{
    void setModelImageBuffer(ndarray::Array<Pixel,1,1> const & buffer);
    void setLinearParameterDerivativeBuffer(
        ndarray::Array<Pixel,2,1> const & buffer
    );
    void setNonlinearParameterDerivativeBuffer(
        ndarray::Array<Pixel,2,1> const & buffer
    );

    //@} 
    
    enum ProductFlag {
        MODEL_IMAGE = 1<<0,
        LINEAR_PARAMETER_DERIVATIVE = 1<<1,
        NONLINEAR_PARAMETER_DERIVATIVE = 1<<2,
    };

    int _validProducts;

    Model::ConstPtr _model;
    boost::shared_ptr<lsst::afw::detection::Footprint const> _footprint;
    lsst::afw::image::Wcs::ConstPtr _wcs;
    // Transform from global coordinates to this projection
    lsst::afw::geom::AffineTransform _skyToPixelTransform; 

    ndarray::Array<Pixel, 1, 1> _modelImage;
    ndarray::Array<Pixel, 2, 1> _linearParameterDerivative;
    ndarray::Array<Pixel, 2, 1> _nonlinearParameterDerivative;
};

}}} // namespace lsst::meas::multifit

#endif // !LSST_MEAS_MULTIFIT_MODEL_PROJECTION_H
