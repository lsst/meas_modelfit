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
 
#include "lsst/meas/multifit/ComponentModelProjection.h"
#include "lsst/meas/multifit/matrices.h"

namespace multifit = lsst::meas::multifit;
namespace afwDet = lsst::afw::detection;
namespace afwImg = lsst::afw::image;
namespace afwGeom = lsst::afw::geom;
namespace afwMath = lsst::afw::math;
/**
 * Construct a ComponentModelProjection
 *
 */
multifit::ComponentModelProjection::ComponentModelProjection(
        ComponentModel::ConstPtr const & model,
        lsst::afw::detection::Psf::ConstPtr const & psf,
        lsst::afw::geom::AffineTransform const & transform,
        CONST_PTR(lsst::afw::detection::Footprint) const & footprint
): ModelProjection(model, transform, footprint),
    _validProducts(0),
    _morphologyProjection(), 
    _translationDerivative(), 
    _projectedParameterDerivative() 
{
    afwMath::Kernel::ConstPtr kernel = psf->getKernel();

    _morphologyProjection = model->getMorphology()->makeProjection(
        afwGeom::ExtentI::make(
            kernel->getWidth(), kernel->getHeight()
        ), 
        getTransform()
    );
}

void multifit::ComponentModelProjection::_computeNonlinearParameterDerivative(
    ndarray::Array<Pixel,2,1> const & matrix
) {
    MatrixMap matrixMap(
        matrix.getData(),
        matrix.getStride<0>(),
        matrix.getSize<0>()
    );
    int nAstrometry = getModel()->getAstrometry()->getParameterSize();
    if (hasTranslationDerivative()) {
        _ensureTranslationDerivative();
        // TODO: Move this into an inline function when possible.
        MatrixMapBlock astrometryView(
            matrixMap, 
            0, 0, 
            matrix.getSize<1>(), nAstrometry
        );
        // END TODO
        TranslationMatrixMap translationView(getTranslationMatrixView());
        components::Astrometry::DerivativeMatrix astrometryDerivative(
            getModel()->getAstrometry()->differentiate()
        );
        astrometryView = translationView * 
            getTransform().getLinear().getMatrix() * astrometryDerivative;
    }
    if (hasProjectedParameterDerivative()) {
        _ensureProjectedParameterDerivative();
        // TODO: Move this into an inline function when possible.
        MatrixMapBlock morphologyView(
            matrixMap, 
            0, nAstrometry, 
            matrix.getSize<1>(), 
            getMorphologyProjection()->getNonlinearParameterSize()
        );


        // END TODO
        // TODO: Move this into an inline function when possible.
        MatrixMap projectedMap(
            _projectedParameterDerivative.getData(),
            _projectedParameterDerivative.getSize<1>(), 
            _projectedParameterDerivative.getSize<0>()
        );
        // END TODO
        morphologyView = projectedMap * 
            (*_morphologyProjection->computeProjectedParameterJacobian());
    }
}

/**
 * Ensure's that _translationDerivative is up to date
 *
 * If _translationDerivative array has not been allocated, do so first.
 * If it is not up to date, call _computeTranslationParameterDerivative
 */
void multifit::ComponentModelProjection::_ensureTranslationDerivative() {
    if (_translationDerivative.empty()) {
        ndarray::shallow(_translationDerivative) = ndarray::allocate<Allocator>(
            ndarray::makeVector(2, getFootprint()->getNpix())
        );
        _validProducts &= (~TRANSLATION_DERIVATIVE);
    }
    if (!(_validProducts & TRANSLATION_DERIVATIVE)) {
        _computeTranslationDerivative(_translationDerivative);
        _validProducts |= TRANSLATION_DERIVATIVE;
    }
}

/**
 * Ensure that _projectedParameterDerivative is up to date
 *
 * If _projectedParameterDerivative array has not been allocated, do so first.
 * If it is not up to date, call _computeProjectedParameterDerivative
 */
void multifit::ComponentModelProjection::_ensureProjectedParameterDerivative() {
    if (_projectedParameterDerivative.empty()) {
        ndarray::shallow(_projectedParameterDerivative) =  ndarray::allocate<Allocator>(
            ndarray::makeVector(
                getMorphologyProjection()->getNonlinearParameterSize(),
                getFootprint()->getNpix()
            )
        );
        _validProducts &= (~PROJECTED_PARAMETER_DERIVATIVE);
    }
    if (!(_validProducts & PROJECTED_PARAMETER_DERIVATIVE)) {
        _computeProjectedParameterDerivative(_projectedParameterDerivative);
        _validProducts |= PROJECTED_PARAMETER_DERIVATIVE;
    }
}
