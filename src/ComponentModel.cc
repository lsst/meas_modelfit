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
 
#include "lsst/meas/multifit/ComponentModel.h"
#include "lsst/meas/multifit/FourierModelProjection.h"
#include "lsst/meas/multifit/footprintUtils.h"

namespace multifit = lsst::meas::multifit;
namespace afwDet = lsst::afw::detection;
namespace afwImg = lsst::afw::image;
namespace afwGeom = lsst::afw::geom;
namespace afwMath = lsst::afw::math;
//-- ComponentModel ------------------------------------------------------------

afwDet::Footprint::Ptr multifit::ComponentModel::computeProjectionFootprint(
    afwDet::Psf::ConstPtr const & psf,
    afwGeom::AffineTransform const & transform
) const {
    afwGeom::ellipses::Ellipse::Ptr ellipse(computeBoundingEllipse());
    ellipse->transform(transform).inPlace();

    if (psf) {
        afwMath::Kernel::ConstPtr kernel = psf->getKernel();
        int kernelSize = std::max(kernel->getWidth(), kernel->getHeight());
        ellipse->grow(kernelSize/2+1);
    }
    return makeFootprint(*ellipse);
}

afwGeom::BoxD multifit::ComponentModel::computeProjectionEnvelope(
    afwDet::Psf::ConstPtr const & psf,
    afwGeom::AffineTransform const & transform
) const {
    afwGeom::ellipses::Ellipse::Ptr ellipse(computeBoundingEllipse());
    ellipse->transform(transform).inPlace();
    if (psf) {
        afwMath::Kernel::ConstPtr kernel = psf->getKernel();
        int kernelSize = std::max(kernel->getWidth(), kernel->getHeight());
        ellipse->grow(kernelSize/2+1);
    }
    return ellipse->computeEnvelope();
}
afwGeom::ellipses::Ellipse::Ptr multifit::ComponentModel::computeBoundingEllipse() const {
    return _morphology->computeBoundingEllipseCore()->makeEllipse(getPosition());
}

void multifit::ComponentModel::_handleLinearParameterChange() {
    _morphology->_handleLinearParameterChange();
}

void multifit::ComponentModel::_handleNonlinearParameterChange() {
    _astrometry->_handleParameterChange();
    _morphology->_handleNonlinearParameterChange();
}

multifit::ModelProjection::Ptr multifit::ComponentModel::makeProjection(
    afwDet::Psf::ConstPtr const & psf,
    afwGeom::AffineTransform const & transform,
    CONST_PTR(afwDet::Footprint) const & footprint
) const {
    ModelProjection::Ptr projection(
        new FourierModelProjection(
            boost::static_pointer_cast<ComponentModel const>(shared_from_this()),
            psf, transform, footprint 
        )
    );
    _registerProjection(projection);
    return projection;
}
/**
 * Initialize the ComponentModel 
 * The number of nonlinear parameters in a ComponentModel
 * is determined by its components. 
 *
 * @param astrometry determines the astrometry component 
 * @param morphology determines the morphology componente
 * @param initializeParameters if true, the models parameter vectors will be
 *      initialized from the values in the input arguments
 */
multifit::ComponentModel::ComponentModel(
    components::Astrometry::ConstPtr const & astrometry,
    components::Morphology::ConstPtr const & morphology,
    bool initializeParameters
) : Model(
        morphology->getLinearParameterSize(), 
        astrometry->getParameterSize() + morphology->getNonlinearParameterSize()
    )     
{
    _initializeFromComponents(astrometry, morphology, initializeParameters);
}

/**
 * Construct a deep copy of ComponentModel. 
 *
 * ModelProjections associated with model will not be associated with the
 * new copy.
 *
 * @sa clone
 */
multifit::ComponentModel::ComponentModel(
    ComponentModel const & model
) : Model(model) {
    _initializeFromComponents(
        model._astrometry, 
        model._morphology, 
        true
    );
}

/**
 * Shared initialization code for ComponentModel constructors.
 *
 * Handles the cloning of the components, feeding them the appropriate iterators
 * into the _nonlinearParameterVector
 */
void multifit::ComponentModel::_initializeFromComponents(
    components::Astrometry::ConstPtr const & astrometry,
    components::Morphology::ConstPtr const & morphology,
    bool initializeParameters
) {
    if(initializeParameters) {
        int nAstrometry = astrometry->getParameterSize();
        std::copy(
            astrometry->begin(), 
            astrometry->end(), 
            _nonlinearParameters->data()
        );
        std::copy(
            morphology->beginNonlinear(), 
            morphology->endNonlinear(), 
            _nonlinearParameters->data()+nAstrometry
        );
        *_linearParameters << (morphology->getLinearParameters());
    }
    
    _astrometry = astrometry->create(_nonlinearParameters);
    _morphology = morphology->create(
        _linearParameters, 
        _nonlinearParameters, 
        _astrometry->getParameterSize()
    );
}
