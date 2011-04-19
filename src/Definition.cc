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

void Object::requireEllipse() const {
    if (!radius || !ellipticity) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            (boost::format("Object %d lacks an ellipticity and/or radius component.") % id).str()
        );
    }
}

} // namespace definition

Definition::Definition(Definition const & other) :
    frames(other.frames), objects(other.objects), wcs(other.wcs) {}

}}} // namespace lsst::meas::multifit
