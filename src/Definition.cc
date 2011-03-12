#include "lsst/meas/multifit/Definition.h"
#include <set>

namespace lsst { namespace meas { namespace multifit {
namespace definition {

Object Object::makeStar(
    ID id, 
    lsst::afw::geom::Point2D const & position, 
    bool isVariable
) {
    Object r(id);
    r.position = boost::make_shared<PositionComponent>(position);
    r.isVariable = isVariable;
    return r;
}

Object Object::makeGalaxy(
    ID id,
    ModelBasis::Ptr const & basis,
    lsst::afw::geom::ellipses::Ellipse const & ellipse
) {
    Object r(id);
    EllipseCore core(ellipse.getCore());
    r.position = boost::make_shared<PositionComponent>(ellipse.getCenter());
    r.ellipticity = boost::make_shared<EllipticityComponent>(core.getEllipticity());
    r.radius = boost::make_shared<RadiusComponent>(core.getRadius());
    r.basis = basis;
    return r;
}

} // namespace definition

Definition::Definition(Definition const & other) :
    frames(other.frames), objects(other.objects), wcs(other.wcs) {}

}}} // namespace lsst::meas::multifit
