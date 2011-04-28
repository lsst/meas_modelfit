#include "lsst/meas/multifit/grid/Grid.h"
#include <boost/format.hpp>
#include <limits>
#include "lsst/ndarray/eigen.h"

namespace lsst { namespace meas { namespace multifit { namespace grid {

static double const EPSILON = std::sqrt(std::numeric_limits<double>::epsilon());

template <ParameterType E>
std::ostream & operator<<(std::ostream & os, ParameterComponent<E> const & component) {
    return os << (*definition::ParameterComponent<E>::make(component.getValue(), component.isActive()))
              << " [" << component.offset << "]";
}

template std::ostream & operator<<(std::ostream &, ParameterComponent<POSITION> const &);
template std::ostream & operator<<(std::ostream &, ParameterComponent<RADIUS> const &);
template std::ostream & operator<<(std::ostream &, ParameterComponent<ELLIPTICITY> const &);

Object::Object(definition::Object const & def, int coefficientOffset, int frameCount, int filterCount) :
    detail::ObjectBase(def), 
    _coefficientOffset(coefficientOffset),
    _coefficientCount(1)
{
    if (getBasis()) {
        _coefficientCount = getBasis()->getSize();
    }
    if (isVariable()) {
        _coefficientCount *= frameCount;
    } else {
        _coefficientCount *= filterCount;
    }
}

void Object::validate() const {
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
                "Objects with a basis must have a radius and ellipticity component."
            );
        }
    } else {
        if (getRadius() || getEllipticity()) {
            throw LSST_EXCEPT(
                lsst::meas::multifit::InvalidDefinitionError,
                "Objects without a basis cannot have a radius or ellipticity component."
            );
        }
    }
}

lsst::afw::geom::Point2D Object::makePoint(double const * paramIter) const {
    lsst::afw::geom::Point2D result = getPosition()->getValue();
    if (getPosition()->isActive()) {
        double const * p = paramIter + getPosition()->offset;
        result += lsst::afw::geom::Extent2D(p[0], p[1]);
    }
    return result;
}

std::pair<int,double> Object::perturbPoint(lsst::afw::geom::Point2D & point, int n) const {
    if (!getPosition()->isActive()) return std::pair<int,double>(-1, 0.0);
    double parameter = point[n] - getPosition()->getValue()[n];
    std::pair<int,double> result(getPosition()->offset + n, ((parameter < 0) ? EPSILON : -EPSILON));
    point[n] += result.second;
    return result;
}

lsst::afw::geom::ellipses::Ellipse Object::makeEllipse(double const * paramIter) const {
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

std::pair<int,double> Object::perturbEllipse(lsst::afw::geom::ellipses::Ellipse & ellipse, int n) const {
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

void Object::unperturbEllipse(
    lsst::afw::geom::ellipses::Ellipse & ellipse, 
    int n, double perturbation
) const {
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

}}}} // namespace lsst::meas::multifit::grid
