#include "lsst/meas/multifit/Definition.h"
#include <set>

namespace lsst { namespace meas { namespace multifit {
namespace definition {

Object Object::makeStar(
    ID id, 
    lsst::afw::geom::Point2D const & position, 
    bool is_variable
) {
    Object r(id);
    r.position = boost::make_shared<PositionComponent>(position);
    r.is_variable = is_variable;
    return r;
}

Object Object::makeGalaxy(
    ID id,
    modeling::EllipseBasis::Ptr const & basis,
    agl::Ellipse const & ellipse
) {
    Object r(id);
    EllipseCore core(ellipse.getCore());
    r.position = boost::make_shared<PositionComponent>(ellipse.getCenter());
    r.ellipticity = boost::make_shared<EllipticityComponent>(core.getEllipticity());
    r.radius = boost::make_shared<RadiusComponent>(core.getRadius());
    r.basis = basis;
    return r;
}

Filter * Filter::get(std::string const & name) {
    typedef definition::Set<Filter,std::string,&Filter::name> Registry;
    static Registry registry;
    Registry::iterator i = registry.find(name);
    if (i == registry.end()) {
        i = registry.insert(registry.end(), Filter(name));
    }
    return &(*i);
}

} // namespace definition

Definition::Definition(Definition const & other) :
    frames(other.frames), objects(other.objects), wcs(other.wcs) {}

}}} // namespace lsst::meas::multifit
