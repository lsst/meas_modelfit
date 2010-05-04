// -*- lsst-c++ -*-
#include "lsst/meas/multifit/ComponentModel.h"
#include "lsst/meas/multifit/FourierModelProjection.h"
#include "lsst/meas/multifit/footprintUtils.h"

namespace multifit = lsst::meas::multifit;

//-- ComponentModel ------------------------------------------------------------

multifit::Footprint::Ptr multifit::ComponentModel::computeProjectionFootprint(
    PsfConstPtr const & psf,
    WcsConstPtr const & wcs
) const {
    lsst::afw::geom::ellipses::Ellipse::Ptr ellipse(computeBoundingEllipse());
    lsst::afw::geom::PointD center(ellipse->getCenter());
    lsst::afw::geom::AffineTransform wcsTransform(wcs->linearizeAt(center));
    
    ellipse->transform(wcsTransform.invert()).inPlace();
    if (psf) {
        int kernelSize = std::max(psf->getWidth(), psf->getHeight());
        ellipse->grow(kernelSize/2+1);
    }
    return makeFootprint(*ellipse);
}

lsst::afw::geom::BoxD multifit::ComponentModel::computeProjectionEnvelope(
    PsfConstPtr const & psf,
    WcsConstPtr const & wcs
) const {
    lsst::afw::geom::ellipses::Ellipse::Ptr ellipse(computeBoundingEllipse());
    lsst::afw::geom::PointD center(ellipse->getCenter());
    lsst::afw::geom::AffineTransform wcsTransform(wcs->linearizeAt(center));
    
    ellipse->transform(wcsTransform.invert()).inPlace();
    if (psf) {
        int kernelSize = std::max(psf->getWidth(), psf->getHeight());
        ellipse->grow(kernelSize/2+1);
    }
    return ellipse->computeEnvelope();
}

lsst::afw::geom::ellipses::Ellipse::Ptr multifit::ComponentModel::computeBoundingEllipse() const {
    return _morphology->computeBoundingEllipseCore()->makeEllipse(
        _astrometry->computePosition()
    );
}

void multifit::ComponentModel::_handleLinearParameterChange() {
    _morphology->_handleLinearParameterChange();
}

void multifit::ComponentModel::_handleNonlinearParameterChange() {
    _astrometry->_handleParameterChange();
    _morphology->_handleNonlinearParameterChange();
}

multifit::ModelProjection::Ptr multifit::ComponentModel::makeProjection(
    PsfConstPtr const & psf,
    WcsConstPtr const & wcs,
    FootprintConstPtr const & footprint
) const {
    ModelProjection::Ptr projection(
        new FourierModelProjection(
            boost::static_pointer_cast<ComponentModel const>(shared_from_this()),
            psf, wcs, footprint 
        )
    );
    _registerProjection(projection);
    return projection;
}

/**
 * Initialize the ComponentModel, initialize (does not set parameters).
 *
 * The number of nonlinear parameters in a ComponentModel
 * is determined by its components. The ComponentModel also clones the
 * astrometryTemplate and morphologyTemplate for future use.
 *
 * @param linearParameterSize used to initialize linear parameter vector
 * @param astrometryTemplate store clone in order to make sense of nonlinear
 *        parameters
 * @param morphologyTemplate store clone in order to make sense of nonlinear
 *        parameters
 */
multifit::ComponentModel::ComponentModel(
    components::Astrometry::ConstPtr const & astrometry,
    components::Morphology::ConstPtr const & morphology,
    bool initializeParameters=false
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
        false
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
    bool copyParameters
) {
    if(copyParameters) {
        std::copy(astrometry->begin(), astrometry->end(), _nonlinearParameters->data());
        std::copy(
            morphology->beginNonlinear(), 
            morphology->endNonlinear(), 
            _nonlinearParameters->data() + astrometry->getParameterSize()
        );
        *_linearParameters << *(morphology->_getLinearParameters());
    }
    
    _astrometry = astrometry->create(_nonlinearParameters);
    _morphology = morphology->create(
        _linearParameters, 
        _nonlinearParameters, 
        _astrometry->getParameterSize()
    );
}
