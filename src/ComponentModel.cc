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
    if (psf) {
        int kernelSize = std::max(psf->getWidth(), psf->getHeight());
        ellipse->grow(kernelSize/2+1);
    }
    lsst::afw::geom::AffineTransform linearApproximation(
        wcs->getAffineTransform(lsst::afw::geom::convertToImage(ellipse->getCenter()))
    );
    ellipse->transform(linearApproximation.invert());

    return makeFootprint(*ellipse);
}

lsst::afw::geom::BoxD multifit::ComponentModel::computeProjectionEnvelope(
    PsfConstPtr const & psf,
    WcsConstPtr const & wcs
) const {
    lsst::afw::geom::ellipses::Ellipse::Ptr ellipse(computeBoundingEllipse());
    if (psf) {
        int kernelSize = std::max(psf->getWidth(), psf->getHeight());
        ellipse->grow(kernelSize/2+1);
    }

    lsst::afw::geom::AffineTransform linearApproximation(
        wcs->getAffineTransform(lsst::afw::geom::convertToImage(ellipse->getCenter()))
    );
    ellipse->transform(linearApproximation.invert());
    return ellipse->computeEnvelope();
}

lsst::afw::geom::ellipses::Ellipse::Ptr multifit::ComponentModel::computeBoundingEllipse() const {
    return lsst::afw::geom::ellipses::Ellipse::Ptr(
        _morphology->computeBoundingEllipseCore()->makeEllipse(_astrometry->apply())
    );
}

void multifit::ComponentModel::_handleLinearParameterChange() {
    _morphology->_handleLinearParameterChange();
}

void multifit::ComponentModel::_handleNonlinearParameterChange() {
    _astrometry->_handleAstrometryParameterChange();
    _morphology->_handleMorphologyParameterChange();
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

multifit::ComponentModel::ComponentModel(
    int linearParameterSize,
    components::Astrometry::ConstPtr const & astrometryTemplate,
    components::Morphology::ConstPtr const & morphologyTemplate
) : Model(
        linearParameterSize, 
        astrometryTemplate->getAstrometryParameterSize() + morphologyTemplate->getMorphologyParameterSize()
    ) 
{
    _construct(astrometryTemplate, morphologyTemplate);
}

multifit::ComponentModel::ComponentModel(ComponentModel const & model) : Model(model) {
    _construct(model._astrometry, model._morphology);
    setLinearParameters(model.getLinearParameterIter());
    setNonlinearParameters(model.getNonlinearParameterIter());
}

void multifit::ComponentModel::_construct(
    components::Astrometry::ConstPtr const & astrometryTemplate,
    components::Morphology::ConstPtr const & morphologyTemplate
) {
    ParameterConstIterator i = getNonlinearParameterIter();
    _astrometry = astrometryTemplate->create(i);
    i += _astrometry->getAstrometryParameterSize();
    _morphology = morphologyTemplate->create(_linearParameterVector, i);
}
