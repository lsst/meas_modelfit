#include "lsst/meas/multifit/Grid.h"
#include <boost/format.hpp>
#include <limits>

namespace lsst { namespace meas { namespace multifit {

namespace {

static double const EPSILON = std::sqrt(std::numeric_limits<double>::epsilon());

template <parameters::Enum E>
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

template <parameters::Enum E, typename Iterator>
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
            grid::ParameterComponent<E> & grid_component
                = static_cast< grid::ParameterComponent<E> & >(*r.first->first);
            r.first->second = grid_component.makeDefinition();
        }
    }
    for (Iterator i = first; i != last; ++i) {
        definition::Object & obj = *i;
        if (!(obj.*member)) continue;
        typename Map::iterator j = unique.find(obj.*member);
        obj.*member = j->second;
    }
}

template <parameters::Enum E, typename Iterator>
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
            grid::ParameterComponent<E> & grid_component
                = static_cast< grid::ParameterComponent<E> & >(*r.first->first);
            r.first->second = grid_component.makeDefinition(parameters);
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
        throw std::invalid_argument(
            boost::str(
                boost::format("Object with ID %d not found.") % id
            )
        );
    }
    return *iter1;
}

template Object const & find(Array<Object> const &, ID const);
template Frame const & find(Array<Frame> const &, ID const);

Object::Object(definition::Object const & definition_, int offset, int frame_count, int filter_count) :
    definition::Object(definition_), 
    coefficient_offset(offset),
    coefficient_count(1)
{
    if (!position) {
        throw definition::InvalidDefinitionError("All objects must have a position component.");
    }
    if (basis) {
        coefficient_count = basis->getSize();
        if (!radius || !ellipticity) {
            throw definition::InvalidDefinitionError(
                "Objects with a basis must have a radius and ellipticity component."
            );
        }
    } else {
        if (radius || ellipticity) {
            throw definition::InvalidDefinitionError(
                "Objects without a basis cannot have a radius or ellipticity component."
            );
        }
    }
 }

agl::PointD Object::makePoint(double const * param_iter) const {
    agl::PointD result;
    if (position->active) {
        double const * p = param_iter + getPosition().offset;
        result = position->getReference() + agl::ExtentD(p[0], p[1]);
    } else {
        result = position->getPosition();
    }
    return result;
}

std::pair<int,double> Object::perturbPoint(agl::PointD & point, int n) const {
    if (!position->active) return std::pair<int,double>(-1, 0.0);
    double parameter = point[n] - position->getReference()[n];
    std::pair<int,double> result(getPosition().offset + n, ((parameter < 0) ? EPSILON : -EPSILON));
    point[n] += result.second;
    return result;
}

agl::Ellipse Object::makeEllipse(double const * param_iter) const {
    typedef agl::ellipses::Separable<EllipticityComponent::Value,RadiusComponent::Value> EllipseCore;
    agl::Ellipse result(
        EllipseCore(
            ellipticity->getValue(),
            radius->getValue()
        ),
        position->getPosition()
    );
    if (position->active) {
        double const * p = param_iter + getPosition().offset;
        result.getCenter() = position->getReference() + agl::ExtentD(p[0], p[1]);
    }
    if (ellipticity->active) {
        double const * p = param_iter + getEllipticity().offset;
        static_cast<EllipseCore&>(result.getCore()).getEllipticity().setE1(p[0]);
        static_cast<EllipseCore&>(result.getCore()).getEllipticity().setE2(p[1]);
    }
    if (radius->active) {
        double const * p = param_iter + getRadius().offset;
        static_cast<EllipseCore&>(result.getCore()).setRadius(p[0]);
    }
    result.getCore().scale(radius_factor);
    return result;
}

std::pair<int,double> Object::perturbEllipse(agl::Ellipse & ellipse, int n) const {
    if (n < 3) {
        typedef agl::ellipses::Separable<EllipticityComponent::Value,RadiusComponent::Value> Core;
        Core & core = static_cast<Core &>(ellipse.getCore());
        core.scale(1.0 / radius_factor);
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
        core.scale(radius_factor);
        return result;
    } else {
        return perturbPoint(ellipse.getCenter(), n-3);
    }
}

void Object::unperturbEllipse(agl::Ellipse & ellipse, int n, double perturbation) const {
    if (n < 3) {
        typedef agl::ellipses::Separable<EllipticityComponent::Value,RadiusComponent::Value> Core;
        Core & core = static_cast<Core &>(ellipse.getCore());
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

Frame::Frame(definition::Frame const & definition_, int offset, int filter_index_, int frame_index_) :
    definition::Frame(definition_), pixel_offset(offset), pixel_count(footprint.getArea()), 
    filter_index(filter_index_), frame_index(frame_index_), extra(0) 
{}

void Frame::applyWeights(ndarray::Array<double,2,1> const & matrix) const {
    if (!weights.empty()) {
        matrix.transpose() *= weights;
    }
}

void Frame::applyWeights(ndarray::Array<double,1,0> const & vector) const {
    if (!weights.empty()) {
        vector.deep() *= weights;
    }
}

Source::Source(
    Frame const & frame_, Object const & object_, 
    boost::shared_ptr<agl::wcs::Projection const> const & wcs
) :
    frame(frame_), object(object_), 
    transform()
{
    if (wcs) {
        if (!frame.wcs) {
            throw definition::InvalidDefinitionError(
                "If the definition WCS is set, all frames must have a WCS."
            );
        }
        transform = wcs->to(*frame.wcs).linearize(object.position->getPosition());
    } else {
        if (frame.wcs) {
            throw definition::InvalidDefinitionError(
                "If the definition WCS is not set, inividual frames may not have a WCS."
            );
        }
    }
    if (object.basis) {
        if (frame.psf) {
            basis = object.basis->convolve(*frame.psf, transform(object.position->getPosition()));
        } else {
            basis = object.basis;
        }
    } else {
        if (!frame.psf) {
            throw definition::InvalidDefinitionError(
                "All objects must have a basis if any frames do not have a PSF."
            );
        }
    }
}

} // namespace grid

template <typename ObjectIterator, typename FrameIterator>
void Grid::_initialize(
    ObjectIterator const & object_begin, ObjectIterator const & object_end,
    FrameIterator const & frame_begin, FrameIterator const & frame_end
) {
    objects._last = objects._first = reinterpret_cast<grid::Object*>(_object_data.get());
    frames._last = frames._first = reinterpret_cast<grid::Frame*>(_frame_data.get());
    sources._last = sources._first = reinterpret_cast<grid::Source*>(_source_data.get());
    try {
        int frame_count = 0;
        for (FrameIterator i = frame_begin; i != frame_end; ++i) {
            std::pair<FilterMap::iterator,bool> r = filters.insert(
                std::make_pair(i->filter, filter_count)
            );
            if (r.second) ++filter_count;
            grid::Frame * new_frame = new (frames._last++) grid::Frame(
                *i, pixel_count, r.first->second, frame_count
            );
            pixel_count += new_frame->pixel_count;
            ++frame_count;
        }
        for (ObjectIterator i = object_begin; i != object_end; ++i) {
            grid::Object * new_object 
                = new (objects._last++) grid::Object(*i, coefficient_count, frames.size(), filter_count);
            coefficient_count += new_object->coefficient_count;
            new_object->sources._first = sources._last;
            for (grid::Frame * j = frames._first; j != frames._last; ++j) {
                new (sources._last++) grid::Source(*j, *new_object, wcs);
            }
            new_object->sources._last = sources._last;
        }
        initializeGridComponents(objects._first, objects._last, positions._ptr_vec, parameter_count, 
                                 &definition::Object::position);
        initializeGridComponents(objects._first, objects._last, radii._ptr_vec, parameter_count,
                                 &definition::Object::radius);
        initializeGridComponents(objects._first, objects._last, ellipticities._ptr_vec, parameter_count,
                                 &definition::Object::ellipticity);
        for (ObjectArray::const_iterator i = objects.begin(); i != objects.end(); ++i) {
            if (i->radius) {
                i->getRadius().associated_ellipticities.insert(i->ellipticity);
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
    filter_count(0),
    coefficient_count(0),
    pixel_count(0),
    parameter_count(0),
    wcs(definition.wcs),
    _object_data(new char[sizeof(grid::Object) * definition.objects.size()]),
    _frame_data(new char[sizeof(grid::Frame) * definition.frames.size()]),
    _source_data(new char[sizeof(grid::Source) * definition.frames.size() * definition.objects.size()])
{
    _initialize(
        definition.objects.begin(), definition.objects.end(),
        definition.frames.begin(), definition.frames.end()
    );
}

Grid::Grid(Grid const & other) :
    filter_count(0),
    coefficient_count(0),
    pixel_count(0),
    parameter_count(0),
    wcs(other.wcs),
    _object_data(new char[sizeof(grid::Object) * other.objects.size()]),
    _frame_data(new char[sizeof(grid::Frame) * other.frames.size()]),
    _source_data(new char[sizeof(grid::Source) * other.sources.size()])
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

void Grid::writeParameters(double * param_iter) const {
    for (PositionArray::const_iterator i = positions.begin(); i != positions.end(); ++i) {
        i->writeParameters(param_iter);
    }
    for (RadiusArray::const_iterator i = radii.begin(); i != radii.end(); ++i) {
        i->writeParameters(param_iter);
    }
    for (EllipticityArray::const_iterator i = ellipticities.begin(); i != ellipticities.end(); ++i) {
        i->writeParameters(param_iter);
    }
}

}}} // namespace lsst::meas::multifit
