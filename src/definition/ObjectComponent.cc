#include "lsst/meas/multifit/definition/ObjectComponent.h"

namespace lsst { namespace meas { namespace multifit {

namespace definition {

ObjectComponent ObjectComponent::makeStar(
    ID const id, 
    lsst::afw::geom::Point2D const & position, 
    bool const isVariable,
    bool const isPositionActive
) {
    ObjectComponent r(id);
    r.getPosition() = PositionElement::make(position, isPositionActive);
    r.getFluxGroup() = FluxGroup::make(id, 1.0, isVariable);
    return r;
}

ObjectComponent ObjectComponent::makeGalaxy(
    ID const id,
    ModelBasis::Ptr const & basis,
    lsst::afw::geom::ellipses::Ellipse const & ellipse,
    bool const isEllipticityActive,
    bool const isRadiusActive,
    bool const isPositionActive
) {
    ObjectComponent r(id);
    EllipseCore core(ellipse.getCore());
    r.getPosition() = PositionElement::make(ellipse.getCenter(), isPositionActive);
    r.getEllipticity() = EllipticityElement::make(core.getEllipticity(), isEllipticityActive);
    r.getRadius() = RadiusElement::make(core.getRadius(), isRadiusActive);
    r.getFluxGroup() = FluxGroup::make(id, basis->getSize(), false);
    r.getBasis() = basis;
    return r;
}

std::ostream & operator<<(std::ostream & os, ObjectComponent const & obj) {
    os << "ObjectComponent " << obj.id << "(@" << (&obj) << "):\n";
    if (obj.getPosition()) os << "    " << (*obj.getPosition()) << "\n";
    if (obj.getRadius()) os << "    " << (*obj.getRadius()) << " x " << obj.getRadiusFactor() << "\n";
    if (obj.getEllipticity()) os << "    " << (*obj.getEllipticity()) << "\n";
    return os;
}

}}}} // namespace lsst::meas::multifit::definition
