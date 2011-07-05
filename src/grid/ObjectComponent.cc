#include "lsst/meas/multifit/grid/Grid.h"
#include <boost/format.hpp>
#include <limits>
#include "lsst/ndarray/eigen.h"

namespace lsst { namespace meas { namespace multifit { namespace grid {

static double const EPSILON = std::sqrt(std::numeric_limits<double>::epsilon());

void ObjectComponent::validate() const {
    if (!getPosition()) {
        throw LSST_EXCEPT(
            lsst::meas::multifit::InvalidDefinitionError,
            "All objects must have a position component."
        );
    }
    if (getBasis()) {
        if (!getRadius() || !getEllipticity()) {
            throw LSST_EXCEPT(
                lsst::meas::multifit::InvalidDefinitionError,
                "ObjectComponents with a basis must have a radius and ellipticity component."
            );
        }
    } else {
        if (getRadius() || getEllipticity()) {
            throw LSST_EXCEPT(
                lsst::meas::multifit::InvalidDefinitionError,
                "ObjectComponents without a basis cannot have a radius or ellipticity component."
            );
        }
    }
}

void ObjectComponent::requireEllipse() const {
    if (!getRadius() || !getEllipticity()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            (boost::format("ObjectComponent %d lacks an ellipticity and/or radius component.") % id).str()
        );
    }
}

lsst::afw::geom::Point2D ObjectComponent::makePoint(
    lsst::ndarray::Array<double const,1,1> const & parameters
) const {
    lsst::afw::geom::Point2D result = getPosition()->getValue();
    if (getPosition()->isActive()) {
        double const * p = parameters.getData() + getPosition()->offset;
        result += lsst::afw::geom::Extent2D(p[0], p[1]);
    }
    return result;
}

void ObjectComponent::readPoint(
    lsst::ndarray::Array<double,1,1> const & parameters,
    lsst::afw::geom::Point2D const & point
) const {
    if (getPosition()->isActive()) {
        double * p = parameters.getData() + getPosition()->offset;
        p[0] = point.getX() - getPosition()->getValue().getX();
        p[1] = point.getY() - getPosition()->getValue().getY();
    }
}

std::pair<int,double> ObjectComponent::perturbPoint(lsst::afw::geom::Point2D & point, int n) const {
    if (!getPosition()->isActive()) return std::pair<int,double>(-1, 0.0);
    double parameter = point[n] - getPosition()->getValue()[n];
    std::pair<int,double> result(getPosition()->offset + n, ((parameter < 0) ? EPSILON : -EPSILON));
    point[n] += result.second;
    return result;
}

lsst::afw::geom::ellipses::Ellipse ObjectComponent::makeEllipse(
    lsst::ndarray::Array<double const,1,1> const & parameters
) const {
    double const * paramIter = parameters.getData();
    requireEllipse();
    lsst::afw::geom::Ellipse result(
        EllipseCore(
            getEllipticity()->getValue(),
            getRadius()->getValue()
        ),
        getPosition()->getValue()
    );
    if (getPosition()->isActive()) {
        double const * p = paramIter + getPosition()->offset;
        result.getCenter() = getPosition()->getValue() + lsst::afw::geom::Extent2D(p[0], p[1]);
    }
    if (getEllipticity()->isActive()) {
        double const * p = paramIter + getEllipticity()->offset;
        static_cast<EllipseCore&>(result.getCore()).getEllipticity().setE1(p[0]);
        static_cast<EllipseCore&>(result.getCore()).getEllipticity().setE2(p[1]);
    }
    if (getRadius()->isActive()) {
        double const * p = paramIter + getRadius()->offset;
        static_cast<EllipseCore&>(result.getCore()).setRadius(p[0]);
    }
    result.getCore().scale(getRadiusFactor());
    return result;
}

void ObjectComponent::readEllipse(
    lsst::ndarray::Array<double,1,1> const & parameters,
    lsst::afw::geom::ellipses::Ellipse const & ellipse
) const {
    double * paramIter = parameters.getData();
    requireEllipse();
    readPoint(parameters, ellipse.getCenter());
    EllipseCore core(ellipse.getCore());
    core.scale(1.0 / getRadiusFactor());
    if (getEllipticity()->isActive()) {
        double * p = paramIter + getEllipticity()->offset;
        p[0] = core.getEllipticity().getE1();
        p[1] = core.getEllipticity().getE2();
    }
    if (getRadius()->isActive()) {
        double * p = paramIter + getRadius()->offset;
        p[0] = core.getRadius();
    }
}

std::pair<int,double> ObjectComponent::perturbEllipse(
    lsst::afw::geom::ellipses::Ellipse & ellipse, int n
) const {
    requireEllipse();
    if (n < 3) {
        EllipseCore & core = static_cast<EllipseCore &>(ellipse.getCore());
        core.scale(1.0 / getRadiusFactor());
        double * parameter;
        std::pair<int,double> result;
        if (n == 2) {
            parameter = &((double&)core.getRadius());
            if (getRadius()->isActive()) {
                result.first = getRadius()->offset;
                result.second = -EPSILON;
            } else {
                result.first = -1;
                result.second = 0.0;
            }
        } else {
            parameter = &core.getEllipticity().getComplex().real() + n;
            if (getEllipticity()->isActive()) {
                result.first = getEllipticity()->offset + n;
                result.second = (*parameter > 0) ? -EPSILON : EPSILON;
            } else {
                result.first = -1;
                result.second = 0.0;
            }
        }
        *parameter += result.second;
        core.scale(getRadiusFactor());
        return result;
    } else {
        return perturbPoint(ellipse.getCenter(), n-3);
    }
}

void ObjectComponent::unperturbEllipse(
    lsst::afw::geom::ellipses::Ellipse & ellipse, 
    int n, double perturbation
) const {
    requireEllipse();
    if (n < 3) {
        EllipseCore & core = static_cast<EllipseCore &>(ellipse.getCore());
        double * parameter;
        std::pair<int,double> result;
        if (n == 2) {
            parameter = &((double&)core.getRadius());
        } else {
            parameter = &core.getEllipticity().getComplex().real() + n;
        }
        *parameter -= perturbation;
    } else {
        unperturbPoint(ellipse.getCenter(), n-3, perturbation);
    }
}

std::ostream & operator<<(std::ostream & os, ObjectComponent const & obj) {
    os << "ObjectComponent " << obj.id << "(@" << (&obj) << "):\n";
    if (obj.getPosition()) os << "    " << (*obj.getPosition()) << "\n";
    if (obj.getRadius()) os << "    " << (*obj.getRadius()) << " x " << obj.getRadiusFactor() << "\n";
    if (obj.getEllipticity()) os << "    " << (*obj.getEllipticity()) << "\n";
    return os;
}

}}}} // namespace lsst::meas::multifit::grid
