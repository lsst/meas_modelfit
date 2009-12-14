#include "lsst/meas/multifit/ComponentModel.h"
#include "lsst/meas/multifit/FourierModelProjection.h"

namespace multifit = lsst::meas::multifit;

//-- ComponentModel ------------------------------------------------------------

multifit::Footprint::Ptr multifit::ComponentModel::computeProjectionFootprint(
    Kernel::ConstPtr const & kernel,
    Wcs::ConstPtr const & wcs,
    double photFactor
) const {
    lsst::afw::geom::ellipses::Ellipse::Ptr ellipse(computeBoundingEllipse());
    ellipse->grow(kernel.getKernelSize()/2+1);

    //TODO: need wcs linearize api
    //ellipse->transform(*wcs->linearize(ellipse->getCenter()));

    //TODO::make footprint from ellipse
    return boost::make_shared<Footprint>();
}

lsst::afw::image::BBox multifit::componentModel::computeProjectionEnvelope(
    Kernel::ConstPtr const & kernel,
    Wcs::ConstPtr const & wcs,
    double photFactor
) const {
    lsst::afw::math::ellipses::Ellipse::Ptr ellipse(computeBoundingEllipse());
    ellipse->grow(kernel.getKernelSize()/2+1);
    //TODO::need api for linearizing a wcs solution
    //ellipse->transform(*wcs->linearize(ellipse->getCenter()));
    return ellipse->computeEnvelope();
}

lsst::afw::math::ellipses::Ellipse::Ptr multifit::ComponentModel::computeBoundingEllipse() const {
    return lsst::afw::math::ellipses::Ellipse::Ptr(
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
    Kernel::ConstPtr const & kernel,
    WCS::ConstPtr const & wcs,
    Footprint::ConstPtr const & footprint,
    double photFactor,
    int activeProducts
) const {
    ModelProjection::Ptr projection(
        new FourierModelProjection(
            boost::static_pointer_cast<ComponentModel const>(shared_from_this()),
            kernel, wcs, footprint, photFactor, activeProducts
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
    i += _astrometry->getParameterSize();
    _morphology = morphologyTemplate->create(_linearParameterVector, i);
}
