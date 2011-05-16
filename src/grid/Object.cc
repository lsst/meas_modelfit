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

void Object::requireEllipse() const {
    if (!getRadius() || !getEllipticity()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::LogicErrorException,
            (boost::format("Object %d lacks an ellipticity and/or radius component.") % id).str()
        );
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

void Object::readPoint(double * paramIter, lsst::afw::geom::Point2D const & point) const {
    if (getPosition()->isActive()) {
        double * p = paramIter + getPosition()->offset;
        p[0] = point.getX() - getPosition()->getValue().getX();
        p[1] = point.getY() - getPosition()->getValue().getY();
    }
}

std::pair<int,double> Object::perturbPoint(lsst::afw::geom::Point2D & point, int n) const {
    if (!getPosition()->isActive()) return std::pair<int,double>(-1, 0.0);
    double parameter = point[n] - getPosition()->getValue()[n];
    std::pair<int,double> result(getPosition()->offset + n, ((parameter < 0) ? EPSILON : -EPSILON));
    point[n] += result.second;
    return result;
}

Eigen::Matrix2d Object::extractPointMatrix(Eigen::MatrixXd const & matrix) const {
    Eigen::Matrix2d r;
    if (getPosition()->isActive()) {
        r = matrix.block<2,2>(getPosition()->offset, getPosition()->offset);
    } else {
        r.setZero();
    }
    return r;
}

void Object::insertPointMatrix(Eigen::MatrixXd & full, Eigen::Matrix2d const & block) const {
    if (getPosition()->isActive()) {
        full.block<2,2>(getPosition()->offset, getPosition()->offset) = block;
    }
}

lsst::afw::geom::ellipses::Ellipse Object::makeEllipse(double const * paramIter) const {
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

void Object::readEllipse(double * paramIter, lsst::afw::geom::ellipses::Ellipse const & ellipse) const {
    requireEllipse();
    readPoint(paramIter, ellipse.getCenter());
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

Eigen::Matrix5d Object::extractEllipseMatrix(Eigen::MatrixXd const & matrix) const {
    requireEllipse();
    Eigen::Matrix5d r = Eigen::Matrix5d::Zero();
    if (getPosition()->isActive() && getEllipticity()->isActive()) {
        r.block<2,2>(3,0) = matrix.block<2,2>(getPosition()->offset, getEllipticity()->offset);
        r.block<2,2>(0,3) = matrix.block<2,2>(getEllipticity()->offset, getPosition()->offset);
    }
    if (getPosition()->isActive() && getRadius()->isActive()) {
        r.block<2,1>(3,2) = matrix.block<2,1>(getPosition()->offset, getRadius()->offset);
        r.block<1,2>(2,3) = matrix.block<1,2>(getRadius()->offset, getPosition()->offset);
    }
    if (getEllipticity()->isActive() && getRadius()->isActive()) {
        r.block<2,1>(0,2) = matrix.block<2,1>(getEllipticity()->offset, getRadius()->offset);
        r.block<1,2>(2,0) = matrix.block<1,2>(getRadius()->offset, getEllipticity()->offset);
    }
    if (getPosition()->isActive()) {
        r.block<2,2>(3,3) = matrix.block<2,2>(getPosition()->offset, getPosition()->offset);
    }
    if (getEllipticity()->isActive()) {
        r.block<2,2>(0,0) = matrix.block<2,2>(getEllipticity()->offset, getEllipticity()->offset);
    }
    if (getRadius()->isActive()) {
        r(2, 2) = matrix(getRadius()->offset, getRadius()->offset);
        r.row(2) *= getRadiusFactor();
        r.col(2) *= getRadiusFactor();
    }
    return r;
}

void Object::insertEllipseMatrix(Eigen::MatrixXd & full, Eigen::Matrix5d const & block) const {
    requireEllipse();
    double rf = getRadiusFactor();
    if (getPosition()->isActive()) {
        full.block<2,2>(getPosition()->offset, getPosition()->offset) = block.block<2,2>(3,3);
    }
    if (getEllipticity()->isActive()) {
        full.block<2,2>(getEllipticity()->offset, getEllipticity()->offset) = block.block<2,2>(0,0);
    }
    if (getRadius()->isActive()) {
        full(getRadius()->offset, getRadius()->offset) = block(2, 2) / (rf * rf);
    }
    if (getPosition()->isActive() && getEllipticity()->isActive()) {
        full.block<2,2>(getPosition()->offset, getEllipticity()->offset) = block.block<2,2>(3,0);
        full.block<2,2>(getEllipticity()->offset, getPosition()->offset) = block.block<2,2>(0,3);
    }
    if (getPosition()->isActive() && getRadius()->isActive()) {
        full.block<2,1>(getPosition()->offset, getRadius()->offset) = block.block<2,1>(3,2) / rf;
        full.block<1,2>(getRadius()->offset, getPosition()->offset) = block.block<1,2>(2,3) / rf;
    }
    if (getEllipticity()->isActive() && getRadius()->isActive()) {
        full.block<2,1>(getEllipticity()->offset, getRadius()->offset) = block.block<2,1>(0,2) / rf;
        full.block<1,2>(getRadius()->offset, getEllipticity()->offset) = block.block<1,2>(2,0) / rf;
    }
}

std::pair<int,double> Object::perturbEllipse(lsst::afw::geom::ellipses::Ellipse & ellipse, int n) const {
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

void Object::unperturbEllipse(
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

std::ostream & operator<<(std::ostream & os, Object const & obj) {
    os << "Object " << obj.id << "(@" << (&obj) << ") = {"
       << (obj.isVariable() ? "variable" : "nonvariable") << ", Rx" << obj.getRadiusFactor() << "}:\n";
    if (obj.getPosition()) os << "    " << (*obj.getPosition()) << "\n";
    if (obj.getRadius()) os << "    " << (*obj.getRadius()) << " x " << obj.getRadiusFactor() << "\n";
    if (obj.getEllipticity()) os << "    " << (*obj.getEllipticity()) << "\n";
    return os;
}

double Source::computeFluxMean(
    lsst::ndarray::Array<double const,1,1> const & parameters,
    lsst::ndarray::Array<double const,1,1> const & coefficients
) const {
    if (object.getBasis()) {
        Ellipse ellipse(object.makeEllipse(parameters.getData()));
        ndarray::Array<Pixel,1,1> integration(ndarray::allocate(object.getBasis()->getSize()));
        object.getBasis()->integrate(integration);
        return (ellipse.getCore().getArea() / M_PI)
            * ndarray::viewAsEigen(
                coefficients[
                    ndarray::view(getCoefficientOffset(), getCoefficientOffset() + getCoefficientCount())
                ]).dot(ndarray::viewAsEigen(integration));
    }
    return coefficients[getCoefficientOffset()];
}

double Source::computeFluxVariance(
    lsst::ndarray::Array<double const,1,1> const & parameters,
    lsst::ndarray::Array<double const,2,1> const & covariance
) const {
    if (object.getBasis()) {
        Eigen::MatrixXd sigma = ndarray::viewAsEigen(covariance).block(
            getCoefficientOffset(),
            getCoefficientOffset(),
            getCoefficientCount(),
            getCoefficientCount()
        );
        Ellipse ellipse(object.makeEllipse(parameters.getData()));
        ndarray::Array<Pixel,1,1> integration(ndarray::allocate(object.getBasis()->getSize()));
        object.getBasis()->integrate(integration);
        double result = ellipse.getCore().getArea() / M_PI; result *= result;
        result *= ndarray::viewAsEigen(integration).dot(sigma * ndarray::viewAsEigen(integration));
        return result;
    }
    return covariance[getCoefficientOffset()][getCoefficientOffset()];
}

}}}} // namespace lsst::meas::multifit::grid
