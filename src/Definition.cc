#include "lsst/meas/multifit/Definition.h"
#include <set>

namespace lsst { namespace meas { namespace multifit {
namespace definition {

Object Object::makeStar(
    ID id, 
    lsst::afw::geom::Point2D const & position, 
    bool isVariable,
    bool isPositionActive
) {
    Object r(id);
    r.position = boost::make_shared<PositionComponent>(position);
    r.position->active = isPositionActive;
    r.isVariable = isVariable;    
    return r;
}

Object Object::makeGalaxy(
    ID id,
    ModelBasis::Ptr const & basis,
    lsst::afw::geom::ellipses::Ellipse const & ellipse,
    bool isEllipticityActive,
    bool isRadiusActive,
    bool isPositionActive
) {
    Object r(id);
    EllipseCore core(ellipse.getCore());
    r.position = boost::make_shared<PositionComponent>(ellipse.getCenter());
    r.position->active = isPositionActive;
    r.ellipticity = boost::make_shared<EllipticityComponent>(core.getEllipticity());
    r.ellipticity->active = isEllipticityActive;
    r.radius = boost::make_shared<RadiusComponent>(core.getRadius());
    r.radius->active = isRadiusActive;
    r.basis = basis;


    return r;
}

std::ostream & operator<<(std::ostream & os, PositionComponent const & component) {
    os << "< Position [" << component.getReference() << " + " << component.getValue() << "] ";
    if (component.active)
        os << "(active)";
    else
        os << "(inactive)";
    return os << " @ " << (&component) << " >";
}

std::ostream & operator<<(std::ostream & os, RadiusComponent const & component) {
    os << "< Radius [" << component.getValue() << "] ";
    if (component.active)
        os << "(active)";
    else
        os << "(inactive)";
    return os << " @ " << (&component) << " >";
}

std::ostream & operator<<(std::ostream & os, EllipticityComponent const & component) {
    os << "< Ellipticity [" << component.getValue().getComplex() << "] ";
    if (component.active)
        os << "(active)";
    else
        os << "(inactive)";
    return os << " @ " << (&component) << " >";
}

std::ostream & operator<<(std::ostream & os, Object const & obj) {
    os << "< Object " << obj.id;
    if (obj.isVariable)
        os << " (variable)";
    else
        os << " (nonvariable)";
    os << " @ " << (&obj) << ">:\n";
    if (obj.position) os << "    " << (*obj.position) << "\n";
    if (obj.radius) os << "    " << (*obj.radius) << " x " << obj.radiusFactor << "\n";
    if (obj.ellipticity) os << "    " << (*obj.ellipticity) << "\n";
    return os;
}

std::ostream & operator<<(std::ostream & os, Frame const & frame) {
    std::string filterName("undefined");
    try {
        filterName = lsst::afw::image::Filter(frame.filterId).getName();
    } catch (lsst::pex::exceptions::NotFoundException &) {}
    os << "< Frame " << frame.id << " (" << filterName << ")";
    return os << " @ " << (&frame) << " >\n";
}

} // namespace definition

Definition::Definition(Definition const & other) :
    frames(other.frames), objects(other.objects), wcs(other.wcs) {}

}}} // namespace lsst::meas::multifit
