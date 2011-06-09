#include "lsst/meas/multifit/grid/Grid.h"
#include <boost/format.hpp>
#include <limits>
#include "lsst/ndarray/eigen.h"

namespace lsst { namespace meas { namespace multifit { namespace grid {

static Pixel const EPSILON = std::sqrt(std::numeric_limits<Pixel>::epsilon());

template <SharedElementType E>
std::ostream & operator<<(std::ostream & os, SharedElement<E> const & element) {
    return os << (*definition::SharedElement<E>::make(
                      element.getValue(), element.getBounds(), element.isActive()
                  ))
              << " [" << element.offset << "]";
}

template std::ostream & operator<<(std::ostream &, SharedElement<POSITION> const &);
template std::ostream & operator<<(std::ostream &, SharedElement<RADIUS> const &);
template std::ostream & operator<<(std::ostream &, SharedElement<ELLIPTICITY> const &);
template std::ostream & operator<<(std::ostream &, SharedElement<FLUX> const &);

ObjectComponent::ObjectComponent(
    definition::ObjectComponent const & def, int coefficientOffset, int frameCount, int filterCount
) :
    detail::ObjectComponentBase(def), 
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

void ObjectComponent::validate() const {
    if (!getPositionElement()) {
        throw LSST_EXCEPT(
            InvalidDefinitionError,
            "All object components must have a position element."
        );
    }
    if (getBasis()) {
        if (!getRadiusElement() || !getEllipticityElement()) {
            throw LSST_EXCEPT(
                InvalidDefinitionError,
                "ObjectComponents with a basis must have radius and ellipticity elements."
            );
        }
    } else {
        if (getRadiusElement() || getEllipticityElement()) {
            throw LSST_EXCEPT(
                InvalidDefinitionError,
                "ObjectComponents without a basis cannot have a radius or ellipticity element."
            );
        }
    }
}

void ObjectComponent::requireEllipse() const {
    if (!getRadiusElement() || !getEllipticityElement()) {
        throw LSST_EXCEPT(
            pex::exceptions::LogicErrorException,
            (boost::format("ObjectComponent %d lacks an ellipticity and/or radius element.") % id).str()
        );
    }
}

afw::geom::Point2D ObjectComponent::makePoint(ndarray::Array<double const> const & parameters) const {
    afw::geom::Point2D result = getPositionElement()->getValue();
    if (getPositionElement()->isActive()) {
        double const * p = parameters.getData() + getPositionElement()->offset;
        result += afw::geom::Extent2D(p[0], p[1]);
    }
    return result;
}

void ObjectComponent::readPoint(
    ndarray::Array<double,1,1> const & parameters,
    afw::geom::Point2D const & point
) const {
    if (getPositionElement()->isActive()) {
        double * p = parameters.getData() + getPositionElement()->offset;
        p[0] = point.getX() - getPositionElement()->getValue().getX();
        p[1] = point.getY() - getPositionElement()->getValue().getY();
    }
}

std::pair<int,double> ObjectComponent::perturbPoint(afw::geom::Point2D & point, int n) const {
    if (!getPositionElement()->isActive()) return std::pair<int,double>(-1, 0.0);
    double parameter = point[n] - getPositionElement()->getValue()[n];
    std::pair<int,double> result(getPositionElement()->offset + n, ((parameter < 0) ? EPSILON : -EPSILON));
    point[n] += result.second;
    return result;
}

Eigen::Matrix2d ObjectComponent::extractPointMatrix(Eigen::MatrixXd const & matrix) const {
    Eigen::Matrix2d r;
    if (getPositionElement()->isActive()) {
        r = matrix.block<2,2>(getPositionElement()->offset, getPositionElement()->offset);
    } else {
        r.setZero();
    }
    return r;
}

void ObjectComponent::insertPointMatrix(Eigen::MatrixXd & full, Eigen::Matrix2d const & block) const {
    if (getPositionElement()->isActive()) {
        full.block<2,2>(getPositionElement()->offset, getPositionElement()->offset) = block;
    }
}

afw::geom::ellipses::Ellipse ObjectComponent::makeEllipse(
    ndarray::Array<double const,1,1> const & parameters
) const {
    double const * paramIter = parameters.getData();
    requireEllipse();
    afw::geom::Ellipse result(
        EllipseCore(
            getEllipticityElement()->getValue(),
            getRadiusElement()->getValue()
        ),
        getPositionElement()->getValue()
    );
    if (getPositionElement()->isActive()) {
        double const * p = paramIter + getPositionElement()->offset;
        result.getCenter() = getPositionElement()->getValue() + afw::geom::Extent2D(p[0], p[1]);
    }
    if (getEllipticityElement()->isActive()) {
        double const * p = paramIter + getEllipticityElement()->offset;
        static_cast<EllipseCore&>(result.getCore()).getEllipticity().setE1(p[0]);
        static_cast<EllipseCore&>(result.getCore()).getEllipticity().setE2(p[1]);
    }
    if (getRadiusElement()->isActive()) {
        double const * p = paramIter + getRadiusElement()->offset;
        static_cast<EllipseCore&>(result.getCore()).setRadius(p[0]);
    }
    result.getCore().scale(getRadiusFactor());
    return result;
}

void ObjectComponent::readEllipse(
    ndarray::Array<double,1,1> const & paramIter, 
    afw::geom::ellipses::Ellipse const & ellipse
) const {
    requireEllipse();
    readPoint(paramIter, ellipse.getCenter());
    EllipseCore core(ellipse.getCore());
    core.scale(1.0 / getRadiusFactor());
    if (getEllipticityElement()->isActive()) {
        double * p = paramIter + getEllipticityElement()->offset;
        p[0] = core.getEllipticity().getE1();
        p[1] = core.getEllipticity().getE2();
    }
    if (getRadiusElement()->isActive()) {
        double * p = paramIter + getRadiusElement()->offset;
        p[0] = core.getRadius();
    }
}

Eigen::Matrix5d ObjectComponent::extractEllipseMatrix(Eigen::MatrixXd const & matrix) const {
    requireEllipse();
    Eigen::Matrix5d r = Eigen::Matrix5d::Zero();
    if (getPositionElement()->isActive() && getEllipticityElement()->isActive()) {
        r.block<2,2>(3,0) = matrix.block<2,2>(getPositionElement()->offset, getEllipticityElement()->offset);
        r.block<2,2>(0,3) = matrix.block<2,2>(getEllipticityElement()->offset, getPositionElement()->offset);
    }
    if (getPositionElement()->isActive() && getRadiusElement()->isActive()) {
        r.block<2,1>(3,2) = matrix.block<2,1>(getPositionElement()->offset, getRadiusElement()->offset);
        r.block<1,2>(2,3) = matrix.block<1,2>(getRadiusElement()->offset, getPositionElement()->offset);
    }
    if (getEllipticityElement()->isActive() && getRadiusElement()->isActive()) {
        r.block<2,1>(0,2) = matrix.block<2,1>(getEllipticityElement()->offset, getRadiusElement()->offset);
        r.block<1,2>(2,0) = matrix.block<1,2>(getRadiusElement()->offset, getEllipticityElement()->offset);
    }
    if (getPositionElement()->isActive()) {
        r.block<2,2>(3,3) = matrix.block<2,2>(getPositionElement()->offset, getPositionElement()->offset);
    }
    if (getEllipticityElement()->isActive()) {
        r.block<2,2>(0,0) = matrix.block<2,2>(getEllipticityElement()->offset, 
                                              getEllipticityElement()->offset);
    }
    if (getRadiusElement()->isActive()) {
        r(2, 2) = matrix(getRadiusElement()->offset, getRadiusElement()->offset);
        r.row(2) *= getRadiusFactor();
        r.col(2) *= getRadiusFactor();
    }
    return r;
}

void ObjectComponent::insertEllipseMatrix(Eigen::MatrixXd & full, Eigen::Matrix5d const & block) const {
    requireEllipse();
    double rf = getRadiusFactor();
    if (getPositionElement()->isActive()) {
        full.block<2,2>(getPositionElement()->offset, getPositionElement()->offset) = block.block<2,2>(3,3);
    }
    if (getEllipticityElement()->isActive()) {
        full.block<2,2>(getEllipticityElement()->offset, getEllipticityElement()->offset)
            = block.block<2,2>(0,0);
    }
    if (getRadiusElement()->isActive()) {
        full(getRadiusElement()->offset, getRadiusElement()->offset) = block(2, 2) / (rf * rf);
    }
    if (getPositionElement()->isActive() && getEllipticityElement()->isActive()) {
        full.block<2,2>(getPositionElement()->offset, getEllipticityElement()->offset) 
            = block.block<2,2>(3,0);
        full.block<2,2>(getEllipticityElement()->offset, getPositionElement()->offset)
            = block.block<2,2>(0,3);
    }
    if (getPositionElement()->isActive() && getRadiusElement()->isActive()) {
        full.block<2,1>(getPositionElement()->offset, getRadiusElement()->offset) 
            = block.block<2,1>(3,2) / rf;
        full.block<1,2>(getRadiusElement()->offset, getPositionElement()->offset) 
            = block.block<1,2>(2,3) / rf;
    }
    if (getEllipticityElement()->isActive() && getRadiusElement()->isActive()) {
        full.block<2,1>(getEllipticityElement()->offset, getRadiusElement()->offset) 
            = block.block<2,1>(0,2) / rf;
        full.block<1,2>(getRadiusElement()->offset, getEllipticityElement()->offset) 
            = block.block<1,2>(2,0) / rf;
    }
}

std::pair<int,double> ObjectComponent::perturbEllipse(afw::geom::ellipses::Ellipse & ellipse, int n) const {
    requireEllipse();
    if (n < 3) {
        EllipseCore & core = static_cast<EllipseCore &>(ellipse.getCore());
        core.scale(1.0 / getRadiusFactor());
        double * parameter;
        std::pair<int,double> result;
        if (n == 2) {
            parameter = &((double&)core.getRadius());
            if (getRadiusElement()->isActive()) {
                result.first = getRadiusElement()->offset;
                result.second = -EPSILON;
            } else {
                result.first = -1;
                result.second = 0.0;
            }
        } else {
            parameter = &core.getEllipticity().getComplex().real() + n;
            if (getEllipticityElement()->isActive()) {
                result.first = getEllipticityElement()->offset + n;
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
    afw::geom::ellipses::Ellipse & ellipse, 
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
    os << "ObjectComponent " << obj.id << "(@" << (&obj) << ") = {"
       << (obj.isVariable() ? "variable" : "nonvariable") << ", Rx" << obj.getRadiusFactor() << "}:\n";
    if (obj.getPositionElement()) os << "    " << (*obj.getPositionElement()) << "\n";
    if (obj.getRadiusElement()) os << "    " << (*obj.getRadiusElement()) << " x " << obj.getRadiusFactor() << "\n";
    if (obj.getEllipticityElement()) os << "    " << (*obj.getEllipticityElement()) << "\n";
    return os;
}

Frame::Frame(definition::Frame const & def, int pixelOffset, int filterIndex, int frameIndex) :
    detail::FrameBase(def, true), _pixelOffset(pixelOffset), _pixelCount(_footprint->getNpix()), 
    _filterIndex(filterIndex), _frameIndex(frameIndex)
{}

void Frame::applyWeights(ndarray::Array<Pixel,2,1> const & matrix) const {
    if (!_weights.empty()) {
        matrix.deep() *= _weights;
    }
}

void Frame::applyWeights(ndarray::Array<Pixel,1,0> const & vector) const {
    if (!_weights.empty()) {
        vector.deep() *= _weights;
    }
}

std::ostream & operator<<(std::ostream & os, Frame const & frame) {
    std::string filterName("undefined");
    try {
        filterName = afw::image::Filter(frame.getFilterId()).getName();
    } catch (pex::exceptions::NotFoundException &) {}
    return os << "Frame " << frame.id << " (@" << (&frame) << ") = {" << filterName 
              << ", " << frame.getFootprint()->getArea() << "pix}\n";
}

SourceComponent::SourceComponent(
    Frame const & frame_, ObjectComponent const & object_, 
    CONST_PTR(afw::image::Wcs) const & wcs
) :
    frame(frame_), object(object_), 
    _transform()
{
    if (wcs) {
        if (!frame.getWcs()) {
            throw LSST_EXCEPT(
                InvalidDefinitionError,
                "If the definition WCS is set, all frames must have a WCS."
            );
        }
	afw::geom::Point2D point = object.getPositionElement()->getValue();
        _transform = frame.getWcs()->linearizeSkyToPixel(point) * wcs->linearizePixelToSky(point);
    } else {
        if (frame.getWcs()) {
            throw LSST_EXCEPT(
                InvalidDefinitionError,
                "If the definition WCS is not set, inividual frames may not have a WCS."
            );
        }
    }
    if (frame.getPsf()) {
        _localPsf = frame.getPsf()->getLocalPsf(_transform(object.getPositionElement()->getValue()));
    }
    if (object.getBasis()) {
        if (frame.getPsf()) {            
            _basis = object.getBasis()->convolve(_localPsf);
        } else {
            _basis = object.getBasis();
        }
    } else {
        if (!frame.getPsf()) {
            throw LSST_EXCEPT(
                InvalidDefinitionError,
                "All objects must have a basis if any frames do not have a PSF."
            );
        }
    }
}

template <typename T>
T const &
find(Array<T> const & array, ID const id) {
    typename Array<T>::const_iterator iter1 = array.begin();
    typename Array<T>::const_iterator iter2;
    int count, step;
    count = array.size();
    while (count > 0) {
        iter2 = iter1;
        step = count / 2;
        iter2 += step;
        if (iter2->id < id) {
            iter1 = ++iter2;
            count -= step + 1;
        } else {
            count = step;
        }
    }
    if (iter1->id != id) {
        throw LSST_EXCEPT(
            pex::exceptions::InvalidParameterException,
            (boost::format("ObjectComponent or Frame with ID %d not found.") % id).str()
        );
    }
    return *iter1;
}

template ObjectComponent const & find(Array<ObjectComponent> const &, ID const);
template Frame const & find(Array<Frame> const &, ID const);

}}}} // namespace lsst::meas::multifit::grid
