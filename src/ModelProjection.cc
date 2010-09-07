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
 
#include "lsst/meas/multifit/ModelProjection.h"
#include "lsst/meas/multifit/matrices.h"
#include "lsst/pex/exceptions/Runtime.h"

namespace multifit = lsst::meas::multifit;


/**
 * \brief Construct a ModelProjection
 *
 * This projection is associated with model, and will be notified
 * whenever the parameters of model are modified.
 *
 * @throw lsst::pex::exceptions::InvalidParameterException footprint is empty
 */
multifit::ModelProjection::ModelProjection(
    Model::ConstPtr const & model,
    lsst::afw::geom::AffineTransform const & transform,
    boost::shared_ptr<lsst::afw::detection::Footprint const> const & footprint
) : _validProducts(0),
    _model(model),
    _footprint(footprint), 
    _wcs(),
    _transform(new lsst::afw::geom::AffineTransform(transform)),
    _modelImage(), 
    _linearParameterDerivative(), 
    _nonlinearParameterDerivative()
{
    if (_footprint->getNpix() <= 0) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Cannot create model projection with empty footprint."
        );
    }
}
multifit::ModelProjection::ModelProjection(
    Model::ConstPtr const & model,
    CONST_PTR(lsst::afw::image::Wcs) const & wcs,
    boost::shared_ptr<lsst::afw::detection::Footprint const> const & footprint
) : _validProducts(0),
    _model(model),
    _footprint(footprint), 
    _wcs(wcs),
    _transform(),
    _modelImage(), 
    _linearParameterDerivative(), 
    _nonlinearParameterDerivative()
{
    if(wcs) {
        lsst::afw::geom::AffineTransform wcsTransform(
            wcs->linearizeAt(model->computePosition())
        ); 
        _transform = boost::make_shared<lsst::afw::geom::AffineTransform>(
            wcsTransform.invert()            
        );
    }

    if (_footprint->getNpix() <= 0) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Cannot create model projection with empty footprint."
        );
    }
}

/**
 * Compute the image the model given the wcs, and psf of this projection
 */
ndarray::Array<multifit::Pixel const,1,1> multifit::ModelProjection::computeModelImage() {
    if (_modelImage.empty()) {
        ndarray::Array<Pixel, 1, 1> buffer(
            ndarray::allocate<Allocator>(
                ndarray::makeVector(_footprint->getNpix())
            )
        );
        setModelImageBuffer(buffer);

        _validProducts &= (~MODEL_IMAGE);
    }
    if (!(_validProducts & MODEL_IMAGE)) {
        _computeModelImage(_modelImage);
        _validProducts |= MODEL_IMAGE;
    }
    return _modelImage;
}

/**
 * Compute the derivative of the model with repect to the linear parameters
 *
 * @throw lsst::pex::exceptions::LogicErrorException if unable to compute a
 * derivative
 */
ndarray::Array<multifit::Pixel const,2,1> multifit::ModelProjection::computeLinearParameterDerivative() {
    if (!hasLinearParameterDerivative()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Cannot compute linear parameter derivative."
        );
    }
    if (_linearParameterDerivative.empty()) {
        ndarray::Array<Pixel, 2, 1> buffer(
            ndarray::allocate<Allocator>(
                ndarray::makeVector(getLinearParameterSize(),_footprint->getNpix())
            )
        );
        setLinearParameterDerivativeBuffer(buffer);

        _validProducts &= (~LINEAR_PARAMETER_DERIVATIVE);
    }
    if (!(_validProducts & LINEAR_PARAMETER_DERIVATIVE)) {
        _computeLinearParameterDerivative(_linearParameterDerivative);
        _validProducts |= LINEAR_PARAMETER_DERIVATIVE;
    }
    return _linearParameterDerivative;
}
/**
 * Compute the derivative of the model with repect to the nonlinear parameters
 *
 * @throw lsst::pex::exceptions::LogicErrorException if unable to compute a
 * derivative
 */
ndarray::Array<multifit::Pixel const,2,1> multifit::ModelProjection::computeNonlinearParameterDerivative() {
    if (!hasNonlinearParameterDerivative()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Cannot compute nonlinear parameter derivative."
        );
    }
    if (_nonlinearParameterDerivative.empty()) {
        ndarray::Array<Pixel, 2, 1> buffer(
            ndarray::allocate<Allocator>(
                ndarray::makeVector(getNonlinearParameterSize(),_footprint->getNpix())
            )
        );
        setNonlinearParameterDerivativeBuffer(buffer);
        _validProducts &= (~NONLINEAR_PARAMETER_DERIVATIVE);
    }
    if (!(_validProducts & NONLINEAR_PARAMETER_DERIVATIVE)) {
        _computeNonlinearParameterDerivative(_nonlinearParameterDerivative);
        _validProducts |= NONLINEAR_PARAMETER_DERIVATIVE;
    }
    return _nonlinearParameterDerivative;
}

/**
 * Set the buffer in which the model image will be computed when
 * computeModelImage is invoked.
 *
 * If computeModelImage is called before this function, then computeModelImage
 * will allocate its own workspace.
 *
 * @throw lsst::pex::exceptions::InvalidParameterException if size of argument
 * buffer does not match number of pixels in footprint
 */
void multifit::ModelProjection::setModelImageBuffer(ndarray::Array<Pixel,1,1> const & buffer) {
    if (buffer.getSize<0>() != _footprint->getNpix()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Model image buffer's size must match Footprint size."
        );
    }
    ndarray::shallow(_modelImage) = buffer;
    _modelImage=0.0;
    _validProducts &= (~MODEL_IMAGE);
}
/**
 * Set the buffer in which the derivative with respect to linear parameters will
 * be computed when computeLinearParameterDerivative is invoked.
 *
 * If computeLinearParameterDerivative is called before this function, then 
 * computeLinearParameterDerivative will allocate its own workspace.
 *
 * @throw lsst::pex::exceptions::LogicErrorException if unable to compute a
 * derivative
 *
 * @throw lsst::pex::exceptions::InvalidParameterException if first dimension
 * of argument buffer does not match number of linear parameters, or if the
 * second dimension does not match the number of pixels in footprint
 */
void multifit::ModelProjection::setLinearParameterDerivativeBuffer(ndarray::Array<Pixel,2,1> const & buffer) {
    if (!hasLinearParameterDerivative()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Projection has no linear parameter derivative."
        );
    }
    if (buffer.getSize<1>() != _footprint->getNpix()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Linear parameter derivative buffer's image dimension must match "
            "Footprint size."
        );
    }
    if (buffer.getSize<0>() != getLinearParameterSize()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Linear parameter derivative buffer's parameter dimension must "
            "match linear parameter size."
        );
    }
    ndarray::shallow(_linearParameterDerivative) = buffer;
    _linearParameterDerivative = 0.0;
    _validProducts &= (~LINEAR_PARAMETER_DERIVATIVE);
}
/**
 * Set the buffer in which the derivative with respec to linear parameters will
 * be computed when computeNonlinearParameterDerivative is invoked.
 *
 * If computeNonlinearParameterDerivative is called before this function, then 
 * computeNonlinearParameterDerivative will allocate its own workspace.
 *
 * @throw lsst::pex::exceptions::LogicErrorException if unable to compute a
 * derivative
 * 
 * @throw lsst::pex::exceptions::InvalidParameterException if first dimension
 * of argument buffer does not match number of nonlinear parameters, or if the
 * second dimension does not match the number of pixels in footprint
 */
void multifit::ModelProjection::setNonlinearParameterDerivativeBuffer(
    ndarray::Array<Pixel,2,1> const & buffer
) {
    if (!hasNonlinearParameterDerivative()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            "Projection has no nonlinear parameter derivative."
        );
    }
    if (buffer.getSize<1>() != _footprint->getNpix()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Nonlinear parameter derivative buffer's image dimension must match "
            "Footprint size."
        );
    }
    if (buffer.getSize<0>() != getNonlinearParameterSize()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            "Nonlinear parameter derivative buffer's parameter dimension must "
            "match nonlinear parameter size."
        );
    }
    ndarray::shallow(_nonlinearParameterDerivative) = buffer;
    _nonlinearParameterDerivative = 0.0;
    _validProducts &= (~NONLINEAR_PARAMETER_DERIVATIVE);
}


void multifit::ModelProjection::_computeModelImage(ndarray::Array<Pixel,1,1> const & vector) {
    ndarray::Array<Pixel const,2,1> array(computeLinearParameterDerivative());
    MatrixMap map(array.getData(), array.getStride<0>(), array.getSize<0>());
    MatrixMapBlock block(map, 0, 0, array.getSize<1>(), array.getSize<0>());
    getVectorView(vector) = block * (_model->getLinearParameters());
}

void multifit::ModelProjection::_handleLinearParameterChange() {
    _validProducts &= (~MODEL_IMAGE);
    _validProducts &= (~NONLINEAR_PARAMETER_DERIVATIVE);
}

void multifit::ModelProjection::_handleNonlinearParameterChange() {
    _validProducts &= (~MODEL_IMAGE);
    _validProducts &= (~LINEAR_PARAMETER_DERIVATIVE);
    _validProducts &= (~NONLINEAR_PARAMETER_DERIVATIVE);
}

