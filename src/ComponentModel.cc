#include "lsst/meas/multifit/ComponentModel.h"
#include "lsst/meas/multifit/FourierModelProjection.h"

namespace multifit = lsst::meas::multifit;

//-- ComponentModel ------------------------------------------------------------

multifit::Footprint::Ptr multifit::ComponentModel::computeProjectionFootprint(
    KernelConstPtr const & kernel,
    WcsConstPtr const & wcs
) const {
    lsst::afw::geom::ellipses::Ellipse::Ptr ellipse(computeBoundingEllipse());
    int kernelSize = std::max(kernel->getWidth(), kernel->getHeight());
    ellipse->grow(kernelSize/2+1);

    //TODO: need wcs linearize api
    //ellipse->transform(*wcs->linearize(ellipse->getCenter()));

    //TODO::make footprint from ellipse
    return boost::make_shared<Footprint>();
}

lsst::afw::geom::Box2D multifit::ComponentModel::computeProjectionEnvelope(
    KernelConstPtr const & kernel,
    WcsConstPtr const & wcs
) const {
    lsst::afw::geom::ellipses::Ellipse::Ptr ellipse(computeBoundingEllipse());
    int kernelSize = std::max(kernel->getWidth(), kernel->getHeight());
    ellipse->grow(kernelSize/2+1);
    //TODO::need api for linearizing a wcs solution
    //ellipse->transform(*wcs->linearize(ellipse->getCenter()));
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
    KernelConstPtr const & kernel,
    WcsConstPtr const & wcs,
    FootprintConstPtr const & footprint
) const {
    ModelProjection::Ptr projection(
        new FourierModelProjection(
            boost::static_pointer_cast<ComponentModel const>(shared_from_this()),
            kernel, wcs, footprint 
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
