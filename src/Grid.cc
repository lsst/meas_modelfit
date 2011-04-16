#include "lsst/meas/multifit/Grid.h"
#include "lsst/meas/multifit/constants.h"
#include <boost/format.hpp>
#include <limits>
#include "lsst/ndarray/eigen.h"

namespace lsst { namespace meas { namespace multifit {

namespace {

static double const EPSILON = std::sqrt(std::numeric_limits<double>::epsilon());

template <ParameterType E>
void initializeGridComponents(
    grid::Object * const first, grid::Object * const last,
    std::vector< boost::shared_ptr< grid::ParameterComponent<E> > > & vec, int & offset,
    boost::shared_ptr< definition::ParameterComponent<E> > definition::Object::*member
) {
    typedef boost::shared_ptr< definition::ParameterComponent<E> > DPtr;
    typedef boost::shared_ptr< grid::ParameterComponent<E> > GPtr;
    typedef std::map<DPtr,GPtr> Map;
    Map unique;
    for (grid::Object * i = first; i != last; ++i) {
        definition::Object & obj = *i;
        if (!(obj.*member)) continue;
        std::pair<DPtr,GPtr> item(obj.*member, GPtr());
        std::pair<typename Map::iterator,bool> r = unique.insert(item);
        if (r.second) {
            vec.push_back(
                r.first->second = boost::make_shared< grid::ParameterComponent<E> >(*r.first->first, offset)
            );
            if (r.first->second->active) {
                offset += definition::ParameterComponent<E>::SIZE;
            }
        }
    }
    for (grid::Object * i = first; i != last; ++i) {
        definition::Object & obj = *i;
        if (!(obj.*member)) continue;
        typename Map::iterator j = unique.find(obj.*member);
        obj.*member = j->second;
    }
}

template <ParameterType E, typename Iterator>
void initializeDefinitionComponents(
    Iterator const & first, Iterator const & last,
    boost::shared_ptr< definition::ParameterComponent<E> > definition::Object::*member
) {
    typedef boost::shared_ptr< definition::ParameterComponent<E> > Ptr;
    typedef std::map<Ptr,Ptr> Map;
    Map unique;
    for (Iterator i = first; i != last; ++i) {
        definition::Object & obj = *i;
        if (!(obj.*member)) continue;
        std::pair<Ptr,Ptr> item(obj.*member, Ptr());
        std::pair<typename Map::iterator,bool> r = unique.insert(item);
        if (r.second) {
            grid::ParameterComponent<E> & gridComponent
                = static_cast< grid::ParameterComponent<E> & >(*r.first->first);
            r.first->second = gridComponent.makeDefinition();
        }
    }
    for (Iterator i = first; i != last; ++i) {
        definition::Object & obj = *i;
        if (!(obj.*member)) continue;
        typename Map::iterator j = unique.find(obj.*member);
        obj.*member = j->second;
    }
}

template <ParameterType E, typename Iterator>
void initializeDefinitionComponents(
    Iterator const & first, Iterator const & last, double const * parameters,
    boost::shared_ptr< definition::ParameterComponent<E> > definition::Object::*member
) {
    typedef boost::shared_ptr< definition::ParameterComponent<E> > Ptr;
    typedef std::map<Ptr,Ptr> Map;
    Map unique;
    for (Iterator i = first; i != last; ++i) {
        definition::Object & obj = *i;
        if (!(obj.*member)) continue;
        std::pair<Ptr,Ptr> item(obj.*member, Ptr());
        std::pair<typename Map::iterator,bool> r = unique.insert(item);
        if (r.second) {
            grid::ParameterComponent<E> & gridComponent
                = static_cast< grid::ParameterComponent<E> & >(*r.first->first);
            r.first->second = gridComponent.makeDefinition(parameters);
        }
    }
    for (Iterator i = first; i != last; ++i) {
        definition::Object & obj = *i;
        if (!(obj.*member)) continue;
        typename Map::iterator j = unique.find(obj.*member);
        obj.*member = j->second;
    }
}

template <typename T>
struct DestroyGridElement {
    void operator()(T & p) const { p.~T(); }
};

} // unnamed

namespace grid {

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
            lsst::pex::exceptions::InvalidParameterException,
            (boost::format("Object with ID %d not found.") % id).str()
        );
    }
    return *iter1;
}

template Object const & find(Array<Object> const &, ID const);
template Frame const & find(Array<Frame> const &, ID const);

Object::Object(definition::Object const & definition_, int offset, int frameCount, int filterCount) :
    definition::Object(definition_), 
    coefficientOffset(offset),
    coefficientCount(1)
{
    if (!position) {
        throw LSST_EXCEPT(
            lsst::meas::multifit::definition::InvalidDefinitionError,
            "All objects must have a position component."
        );
    }
    if (basis) {
        coefficientCount = basis->getSize();
        if (!radius || !ellipticity) {
            throw LSST_EXCEPT(
                lsst::meas::multifit::definition::InvalidDefinitionError,
                "Objects with a basis must have a radius and ellipticity component."
            );
        }
    } else {
        if (radius || ellipticity) {
            throw LSST_EXCEPT(
                lsst::meas::multifit::definition::InvalidDefinitionError,
                "Objects without a basis cannot have a radius or ellipticity component."
            );
        }
    }
 }

lsst::afw::geom::Point2D Object::makePoint(double const * paramIter) const {
    lsst::afw::geom::Point2D result;
    if (position->active) {
        double const * p = paramIter + getPosition().offset;
        result = position->getReference() + lsst::afw::geom::Extent2D(p[0], p[1]);
    } else {
        result = position->getPosition();
    }
    return result;
}

std::pair<int,double> Object::perturbPoint(lsst::afw::geom::Point2D & point, int n) const {
    if (!position->active) return std::pair<int,double>(-1, 0.0);
    double parameter = point[n] - position->getReference()[n];
    std::pair<int,double> result(getPosition().offset + n, ((parameter < 0) ? EPSILON : -EPSILON));
    point[n] += result.second;
    return result;
}

lsst::afw::geom::ellipses::Ellipse Object::makeEllipse(double const * paramIter) const {
    lsst::afw::geom::Ellipse result(
        EllipseCore(
            ellipticity->getValue(),
            radius->getValue()
        ),
        position->getPosition()
    );
    if (position->active) {
        double const * p = paramIter + getPosition().offset;
        result.getCenter() = position->getReference() + lsst::afw::geom::Extent2D(p[0], p[1]);
    }
    if (ellipticity->active) {
        double const * p = paramIter + getEllipticity().offset;
        static_cast<EllipseCore&>(result.getCore()).getEllipticity().setE1(p[0]);
        static_cast<EllipseCore&>(result.getCore()).getEllipticity().setE2(p[1]);
    }
    if (radius->active) {
        double const * p = paramIter + getRadius().offset;
        static_cast<EllipseCore&>(result.getCore()).setRadius(p[0]);
    }
    result.getCore().scale(radiusFactor);
    return result;
}

std::pair<int,double> Object::perturbEllipse(lsst::afw::geom::ellipses::Ellipse & ellipse, int n) const {
    if (n < 3) {
        EllipseCore & core = static_cast<EllipseCore &>(ellipse.getCore());
        core.scale(1.0 / radiusFactor);
        double * parameter;
        std::pair<int,double> result;
        if (n == 2) {
            parameter = &((double&)core.getRadius());
            if (radius->active) {
                result.first = getRadius().offset;
                result.second = -EPSILON;
            } else {
                result.first = -1;
                result.second = 0.0;
            }
        } else {
            parameter = &core.getEllipticity().getComplex().real() + n;
            if (ellipticity->active) {
                result.first = getEllipticity().offset + n;
                result.second = (*parameter > 0) ? -EPSILON : EPSILON;
            } else {
                result.first = -1;
                result.second = 0.0;
            }
        }
        *parameter += result.second;
        core.scale(radiusFactor);
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

Frame::Frame(definition::Frame const & definition_, int offset, int filterIndex_, int frameIndex_) :
    definition::Frame(definition_), pixelOffset(offset), pixelCount(footprint->getNpix()), 
    filterIndex(filterIndex_), frameIndex(frameIndex_), extra(0) 
{}

void Frame::applyWeights(ndarray::Array<double,2,1> const & matrix) const {
    if (!weights.empty()) {
        matrix.deep() *= weights;
    }
}

void Frame::applyWeights(ndarray::Array<double,1,0> const & vector) const {
    if (!weights.empty()) {
        vector.deep() *= weights;
    }
}

Source::Source(
    Frame const & frame_, Object const & object_, 
    CONST_PTR(afw::image::Wcs) const & wcs
) :
    frame(frame_), object(object_), 
    transform()
{
    if(frame.psf) {
        localPsf = frame.psf->getLocalPsf(transform(object.position->getPosition()));
    }
    if (wcs) {
        if (!frame.wcs) {
            throw LSST_EXCEPT(
                lsst::meas::multifit::definition::InvalidDefinitionError,
                "If the definition WCS is set, all frames must have a WCS."
            );
        }
	afw::geom::Point2D point = object.position->getPosition();
        transform = frame.wcs->linearizeSkyToPixel(point)*wcs->linearizePixelToSky(point);
    } else {
        if (frame.wcs) {
            throw LSST_EXCEPT(
                lsst::meas::multifit::definition::InvalidDefinitionError,
                "If the definition WCS is not set, inividual frames may not have a WCS."
            );
        }
    }
    if (object.basis) {
        if (frame.psf) {            
            basis = object.basis->convolve(localPsf);
        } else {
            basis = object.basis;
        }
    } else {
        if (!frame.psf) {
            throw LSST_EXCEPT(
                lsst::meas::multifit::definition::InvalidDefinitionError,
                "All objects must have a basis if any frames do not have a PSF."
            );
        }
    }
}

} // namespace grid

template <typename ObjectIterator, typename FrameIterator>
void Grid::_initialize(
    ObjectIterator const & objectBegin, ObjectIterator const & objectEnd,
    FrameIterator const & frameBegin, FrameIterator const & frameEnd
) {
    objects._last = objects._first = reinterpret_cast<grid::Object*>(_objectData.get());
    frames._last = frames._first = reinterpret_cast<grid::Frame*>(_frameData.get());
    sources._last = sources._first = reinterpret_cast<grid::Source*>(_sourceData.get());
    try {
        int frameCount = 0;
        for (FrameIterator i = frameBegin; i != frameEnd; ++i) {
            std::pair<FilterMap::iterator,bool> r = filters.insert(
                std::make_pair(i->filterId, filterCount)
            );
            if (r.second) ++filterCount;
            grid::Frame * new_frame = new (frames._last++) grid::Frame(
                *i, pixelCount, r.first->second, frameCount
            );
            pixelCount += new_frame->pixelCount;
            ++frameCount;
        }
        for (ObjectIterator i = objectBegin; i != objectEnd; ++i) {
            grid::Object * new_object 
                = new (objects._last++) grid::Object(*i, coefficientCount, frames.size(), filterCount);
            coefficientCount += new_object->coefficientCount;
            new_object->sources._first = sources._last;
            for (grid::Frame * j = frames._first; j != frames._last; ++j) {
                new (sources._last++) grid::Source(*j, *new_object, wcs);
            }
            new_object->sources._last = sources._last;
        }
        initializeGridComponents(objects._first, objects._last, positions._ptrVec, parameterCount, 
                                 &definition::Object::position);
        initializeGridComponents(objects._first, objects._last, radii._ptrVec, parameterCount,
                                 &definition::Object::radius);
        initializeGridComponents(objects._first, objects._last, ellipticities._ptrVec, parameterCount,
                                 &definition::Object::ellipticity);
        for (ObjectArray::const_iterator i = objects.begin(); i != objects.end(); ++i) {
            if (i->radius) {
                i->getRadius().associatedEllipticities.insert(i->ellipticity);
            }
        }
    } catch (...) {
        _destroy();
        throw;
    }
}

void Grid::_destroy() {
    std::for_each(sources._first, sources._last, DestroyGridElement<grid::Source>());
    std::for_each(objects._first, objects._last, DestroyGridElement<grid::Object>());
    std::for_each(frames._first, frames._last, DestroyGridElement<grid::Frame>());
}

Grid::Grid(Definition const & definition) :
    filterCount(0),
    coefficientCount(0),
    pixelCount(0),
    parameterCount(0),
    wcs(definition.wcs),
    _objectData(new char[sizeof(grid::Object) * definition.objects.size()]),
    _frameData(new char[sizeof(grid::Frame) * definition.frames.size()]),
    _sourceData(new char[sizeof(grid::Source) * definition.frames.size() * definition.objects.size()])
{
    _initialize(
        definition.objects.begin(), definition.objects.end(),
        definition.frames.begin(), definition.frames.end()
    );
}

Grid::Grid(Grid const & other) :
    filterCount(0),
    coefficientCount(0),
    pixelCount(0),
    parameterCount(0),
    wcs(other.wcs),
    _objectData(new char[sizeof(grid::Object) * other.objects.size()]),
    _frameData(new char[sizeof(grid::Frame) * other.frames.size()]),
    _sourceData(new char[sizeof(grid::Source) * other.sources.size()])
{
    _initialize(
        other.objects.begin(), other.objects.end(),
        other.frames.begin(), other.frames.end()
    );
}

Definition Grid::makeDefinition(double const * parameters) const {
    Definition r(wcs);
    r.frames.insert(frames.begin(), frames.end());
    r.objects.insert(objects.begin(), objects.end());
    initializeDefinitionComponents(r.objects.begin(), r.objects.end(),
                                   parameters, &definition::Object::position);
    initializeDefinitionComponents(r.objects.begin(), r.objects.end(),
                                   parameters, &definition::Object::radius);
    initializeDefinitionComponents(r.objects.begin(), r.objects.end(),
                                   parameters, &definition::Object::ellipticity);
    return r;
}

Definition Grid::makeDefinition() const {
    Definition r(wcs);
    r.frames.insert(frames.begin(), frames.end());
    r.objects.insert(objects.begin(), objects.end());
    initializeDefinitionComponents(r.objects.begin(), r.objects.end(), &definition::Object::position);
    initializeDefinitionComponents(r.objects.begin(), r.objects.end(), &definition::Object::radius);
    initializeDefinitionComponents(r.objects.begin(), r.objects.end(), &definition::Object::ellipticity);
    return r;
}

void Grid::writeParameters(double * paramIter) const {
    for (PositionArray::const_iterator i = positions.begin(); i != positions.end(); ++i) {
        if (i->active) {
            i->writeParameters(paramIter);
        }
    }
    for (RadiusArray::const_iterator i = radii.begin(); i != radii.end(); ++i) {
        if (i->active) {
            i->writeParameters(paramIter);
        }
    }
    for (EllipticityArray::const_iterator i = ellipticities.begin(); i != ellipticities.end(); ++i) {
        if (i->active) {
            i->writeParameters(paramIter);
        }
    }
}

double Grid::sumLogWeights() const {
    double r = 0;
    for (FrameArray::const_iterator i = frames.begin(); i != frames.end(); ++i) {
        r += ndarray::viewAsEigen(i->weights).cwise().log().sum();
    }
    return r;
}

}}} // namespace lsst::meas::multifit
