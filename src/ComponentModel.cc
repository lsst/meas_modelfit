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
    return lsst::afw::geom::ellipses::Ellipse::Ptr(
        _morphology->computeBoundingEllipseCore()->makeEllipse(_astrometry->computePosition())
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
 * Initialize the ComponentModel (does not set parameters).
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
    int linearParameterSize,
    components::Astrometry::ConstPtr const & astrometryTemplate,
    components::Morphology::ConstPtr const & morphologyTemplate
) : Model(
        linearParameterSize, 
        astrometryTemplate->getParameterSize() + morphologyTemplate->getNonlinearParameterSize()
    ) 
{
    _initializeComponents(astrometryTemplate, morphologyTemplate);
}

/**
 * Construct a deep copy of ComponentModel. 
 *
 * ModelProjections associated with model will not be associated with the
 * new copy.
 *
 * @sa clone
 */
multifit::ComponentModel::ComponentModel(ComponentModel const & model) : Model(model) {
    _initializeComponents(model._astrometry, model._morphology);
    setLinearParameters(model.getLinearParameterIter());
    setNonlinearParameters(model.getNonlinearParameterIter());
}

/**
 * Shared initialization code for ComponentModel constructors.
 *
 * Handles the cloning of the components, feeding them the appropriate iterators
 * into the _nonlinearParameterVector
 */
void multifit::ComponentModel::_initializeComponents(
    components::Astrometry::ConstPtr const & astrometryTemplate,
    components::Morphology::ConstPtr const & morphologyTemplate
) {
    _astrometry = astrometryTemplate->create(_getAstrometryParameterIter());
    _morphology = morphologyTemplate->create(
        _linearParameters,
        _getMorphologyParameterIter()
    );
}
