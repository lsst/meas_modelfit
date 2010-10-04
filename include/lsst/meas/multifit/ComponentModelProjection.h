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
 * Support for viewing and evaluating a ComponentModel
 */
#ifndef LSST_MEAS_MULTIFIT_COMPONENT_MODEL_PROJECTION_H
#define LSST_MEAS_MULTIFIT_COMPONENT_MODEL_PROJECTION_H

#include "lsst/meas/multifit/ComponentModel.h"
#include "lsst/meas/multifit/components/MorphologyProjection.h"
#include "lsst/meas/multifit/matrices.h"

#include "lsst/afw/geom/AffineTransform.h"

namespace lsst {
namespace meas {
namespace multifit {

/**
 *  A projection of a ComponentModel to a particular set of observing 
 *  conditions.
 *  
 *  @sa ModelProjection
 *  @sa ComponentModel
 */
class ComponentModelProjection : public ModelProjection  {
public:
    typedef boost::shared_ptr<ModelProjection> Ptr;
    typedef boost::shared_ptr<ModelProjection const> ConstPtr;

    //static const int WCS_PARAMETER_SIZE = 6;


    /**
     * Immutable access to the ComponentModel this is a projection of
     */
    ComponentModel::ConstPtr getModel() const {
        return boost::static_pointer_cast<ComponentModel const>(
            ModelProjection::getModel()
        );
    }

    /**
     * Immutable reference to the astrometry object this projection is based on.
     */
    CONST_PTR(lsst::meas::multifit::components::Astrometry) getAstrometry() const { 
        return getModel()->getAstrometry(); 
    }

    /**
     * Imutable reference to the MorphologyProjection object this projection is based on.
     */
    lsst::meas::multifit::components::MorphologyProjection::ConstPtr getMorphologyProjection() const {
        return _morphologyProjection;
    }

protected:

    enum ComponentModelProductFlags {
        TRANSLATION_DERIVATIVE = 1 << 0,
        PROJECTED_PARAMETER_DERIVATIVE = 1 << 1
    };
    
    ComponentModelProjection(
        ComponentModel::ConstPtr const & model,
        lsst::afw::detection::Psf::ConstPtr const & psf,
        lsst::afw::geom::AffineTransform const & pixelToPixel,
        CONST_PTR(lsst::afw::detection::Footprint) const & footprint
    );

    /**
     * Compute the image-space (xy) coordinates where the PSF should be centered
     */
    virtual lsst::afw::geom::Point2D _getPsfPosition() const { 
        return (getTransform())(getAstrometry()->computePosition()); 
    }

    /**
     *  @name Protected Product Computers
     *
     *  These are fully implemented by ComponentModel by delegating to _computeTranslationDerivative()
     *  and _computeProjectedParameterDerivative(), and should not generally by reimplemented
     *  by subclasses. Instead, subclasses should implement
     *  _computeTranslationDerivative and _computeProjectedParameterDerivative
     *  as needed
     */
    //@{
    virtual void _computeNonlinearParameterDerivative(ndarray::Array<Pixel,2,1> const & matrix);
#if 0
    //wcs derivatives not yet part of model fitting API
    virtual void _computeWcsParameterDerivative(ndarray::Array<Pixel,2,1> const & matrix);
#endif

    /**
     *  Compute the derivative of the model image with respect to image-coordinate translations
     *  (Footprint-compressed).
     */
    virtual void _computeTranslationDerivative(ndarray::Array<Pixel,2,2> const & matrix) = 0;

    /**
     *  Compute the derivative of the model image with respect to its projected morphology
     *  parameters (Footprint-compressed).
     *
     *  The "projected morphology parameters" correspond to image-frame versions of the
     *  nonlinear morphology parameters, and must transform according to the
     *  MorphologyProjection's getParameterJacobian() and getTransformDerivative() outputs.
     */
    virtual void _computeProjectedParameterDerivative(ndarray::Array<Pixel,2,2> const & matrix) = 0;
    //@}

    /**
     *  @name Product Enabled Checkers
     */
    //@{
    virtual bool hasTranslationDerivative() const { return true; }
    virtual bool hasProjectedParameterDerivative() const {
        return getMorphologyProjection()->getMorphology()->getNonlinearParameterSize() > 0;
    }
    //@}
    
    virtual void _handleLinearParameterChange() {
        ModelProjection::_handleLinearParameterChange();
        _morphologyProjection->_handleLinearParameterChange();
        _validProducts &= ~PROJECTED_PARAMETER_DERIVATIVE;
        _validProducts &= ~TRANSLATION_DERIVATIVE;
    }
    virtual void _handleNonlinearParameterChange() {
        ModelProjection::_handleNonlinearParameterChange();
        _morphologyProjection->_handleNonlinearParameterChange();
        _validProducts &= ~PROJECTED_PARAMETER_DERIVATIVE;
        _validProducts &= ~TRANSLATION_DERIVATIVE;
    }

    /**
     * Mutable reference to the MorphologyProjection object this projection is based on.
     */
    lsst::meas::multifit::components::MorphologyProjection::Ptr _getMorphologyProjection() { 
        return _morphologyProjection;
    }

private:
    typedef Eigen::Map< Eigen::Matrix<Pixel,Eigen::Dynamic,2> > TranslationMatrixMap;

    void _ensureTranslationDerivative();
    void _ensureProjectedParameterDerivative();

    TranslationMatrixMap getTranslationMatrixView() {
        return TranslationMatrixMap(
            _translationDerivative.getData(),
            _translationDerivative.getSize<1>(),
            2
        );
    }

    // enable flags for the products this projection can compute
    int _validProducts;

    // MorphologyProjection this ComponentModelProjection is based on
    lsst::meas::multifit::components::MorphologyProjection::Ptr _morphologyProjection;

    ndarray::Array<Pixel,2,2> _translationDerivative;
    ndarray::Array<Pixel,2,2> _projectedParameterDerivative;
};

}}} // namespace multifit

#endif // !LSST_MEAS_MULTIFIT_COMPONENT_MODEL_PROJECTION_H
