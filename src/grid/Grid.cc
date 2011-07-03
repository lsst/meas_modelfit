#include "lsst/meas/multifit/grid/Grid.h"
#include <boost/format.hpp>

namespace lsst { namespace meas { namespace multifit { namespace grid {

class Initializer {
public:

    typedef std::vector< PTR(ObjectComponent) > ObjectComponentVector;
    typedef std::vector< PTR(SourceComponent) > SourceComponentVector;
    typedef std::vector< PTR(Frame) > FrameVector;
    typedef std::vector< PTR(PositionElement) > PositionVector;
    typedef std::vector< PTR(RadiusElement) > RadiusVector;
    typedef std::vector< PTR(EllipticityElement) > EllipticityVector;
    typedef std::vector< PTR(FluxGroup) > FluxGroupVector;

    int _filterCount;
    int _coefficientCount;
    int _pixelCount;
    int _parameterCount;
    int _constraintCount;

    ObjectComponentVector _objects;
    SourceComponentVector _sources;
    FrameVector _frames;
    FluxGroupVector _groups;
    PositionVector _positions;
    RadiusVector _radii;
    EllipticityVector _ellipticities;

    Grid::FilterMap _filters;

    CONST_PTR(Wcs) _wcs;

    explicit Initializer(Definition const & definition);

    static Definition makeDefinition(Grid const & grid, double const * paramIter);

    template <SharedElementType E>
    void makeGridElements(
        Definition const & def,
        std::map< PTR(definition::SharedElement<E>), PTR(grid::SharedElement<E>) > & map,
        std::vector< PTR(grid::SharedElement<E>) > & vector
    );

    template <SharedElementType E>
    static void transferElementsToDefinition(
        Grid const & input, Definition & output, double const * paramIter
    );
};


Initializer::Initializer(Definition const & def) :
    _filterCount(0), _coefficientCount(0), _pixelCount(0), _parameterCount(0), _constraintCount(0), _wcs()
{
    if (def.getWcs()) _wcs = def.getWcs()->clone();
    _objects.reserve(def.objects.size());
    _frames.reserve(def.frames.size());
    _sources.reserve(def.objects.size() * def.frames.size());

    std::map<definition::PositionElement::Ptr,grid::PositionElement::Ptr> positionMap;
    std::map<definition::RadiusElement::Ptr,grid::RadiusElement::Ptr> radiusMap;
    std::map<definition::EllipticityElement::Ptr,grid::EllipticityElement::Ptr> ellipticityMap;

    makeGridElements(def, positionMap, _positions);
    makeGridElements(def, radiusMap, _radii);
    makeGridElements(def, ellipticityMap, _ellipticities);

    for (Definition::FrameSet::iterator k = def.frames.begin(); k != def.frames.end(); ++k) {
        std::pair<Grid::FilterMap::iterator,bool> r = _filters.insert(
            std::make_pair(k->getFilterId(), _filterCount)
        );
        if (r.second) ++_filterCount;
        _frames.push_back(
            boost::shared_ptr<Frame>(new Frame(*k, _pixelCount, r.first->second, _frames.size()))
        );
        _pixelCount += _frames.back()->getPixelCount();
    }

    std::set<definition::FluxGroup::Ptr> groupsFinished;
    for (Definition::ObjectComponentSet::iterator i = def.objects.begin(); i != def.objects.end(); ++i) {
        if (!i->getFluxGroup()) {
            throw LSST_EXCEPT(
                lsst::meas::multifit::InvalidDefinitionError, 
                (boost::format("ObjectComponent %lld does not have a FluxGroup.") % i->id).str()
            );
        }
        std::pair< std::set<definition::FluxGroup::Ptr>::iterator, bool> r =
            groupsFinished.insert(i->getFluxGroup());
        if (!r.second) continue;
        _groups.push_back(FluxGroup::Ptr(new FluxGroup(*i->getFluxGroup(), _coefficientCount)));
        _groups.back()->components._first = FluxGroup::ComponentArray::iterator(_objects.end());
        for (Definition::ObjectComponentSet::iterator j = i; j != def.objects.end(); ++j) {
            if (i->getFluxGroup() != j->getFluxGroup()) {
                if (i->getFluxGroup()->id == j->getFluxGroup()->id) {
                    throw LSST_EXCEPT(
                        lsst::meas::multifit::InvalidDefinitionError, 
                        (boost::format("Multiple FluxGroups with identical ID %lld.") % i->id).str()
                    );
                } else {
                    continue;
                }
            }
            _objects.push_back(
                boost::shared_ptr<ObjectComponent>(
                    new ObjectComponent(*j, _groups.back()->getSourceCoefficientCount())
                )
            );
            if (j->getPosition())
                _objects.back()->_position = positionMap.find(j->getPosition())->second;
            if (j->getRadius())
                _objects.back()->_radius = radiusMap.find(j->getRadius())->second;
            if (j->getEllipticity())
                _objects.back()->_ellipticity = ellipticityMap.find(j->getEllipticity())->second;
            _groups.back()->_sourceCoefficientCount += _objects.back()->getSourceCoefficientCount();
            _objects.back()->sources._first = ObjectComponent::SourceArray::iterator(_sources.end());
            for (FrameVector::iterator k = _frames.begin(); k != _frames.end(); ++k) {
                _sources.push_back(
                    boost::shared_ptr<SourceComponent>(new SourceComponent(**k, *_objects.back(), _wcs))
                );
            }
            _objects.back()->sources._last = ObjectComponent::SourceArray::iterator(_sources.end());
            _objects.back()->_fluxGroup = _groups.back();
            _objects.back()->validate();
        }
        _groups.back()->components._last = FluxGroup::ComponentArray::iterator(_objects.end());
        _groups.back()->initialize();
        int nFluxInstances = ((_groups.back()->isVariable()) ? _frames.size() : _filterCount);
        _coefficientCount += nFluxInstances * _groups.back()->getSourceCoefficientCount();
        _constraintCount += nFluxInstances * _groups.back()->getConstraintCount();
    }
}

Definition Initializer::makeDefinition(Grid const & grid, double const * paramIter) {
    Wcs::Ptr wcs;
    if (grid.getWcs()) wcs = grid.getWcs()->clone();
    Definition def(wcs);
    for (Grid::FrameArray::iterator i = grid.frames.begin(); i != grid.frames.end(); ++i) {
        def.frames.insert(definition::Frame(*i));
    }
    for (Grid::ObjectComponentArray::iterator i = grid.objects.begin(); i != grid.objects.end(); ++i) {
        def.objects.insert(definition::ObjectComponent(*i));
    }
    transferElementsToDefinition<POSITION>(grid, def, paramIter);
    transferElementsToDefinition<RADIUS>(grid, def, paramIter);
    transferElementsToDefinition<ELLIPTICITY>(grid, def, paramIter);
    return def;
}

template <SharedElementType E>
void Initializer::makeGridElements(
    Definition const & def,
    std::map< PTR(definition::SharedElement<E>), PTR(grid::SharedElement<E>) > & map,
    std::vector< PTR(grid::SharedElement<E>) > & vector
) {
    typedef std::map< PTR(definition::SharedElement<E>), PTR(grid::SharedElement<E>) > Map;
    for (Definition::ObjectComponentSet::iterator i = def.objects.begin(); i != def.objects.end(); ++i) {
        if (!i->template getElement<E>()) continue;
        std::pair<typename Map::iterator,bool> r = map.insert(
            std::make_pair(i->template getElement<E>(), boost::shared_ptr< grid::SharedElement<E> >())
        );
        if (!r.second) continue;
        if (r.first->first->isActive()) {
            r.first->second = boost::shared_ptr< grid::SharedElement<E> >(
                new grid::SharedElement<E>(*r.first->first, _parameterCount)
            );
            _parameterCount += grid::SharedElement<E>::SIZE;
            vector.push_back(r.first->second);
        } else {
            r.first->second = boost::shared_ptr< grid::SharedElement<E> >(
                new grid::SharedElement<E>(*r.first->first, -1)
            );
        }
    }
}

template <SharedElementType E>
void Initializer::transferElementsToDefinition(
    Grid const & input, Definition & output, double const * paramIter
) {
    typedef boost::shared_ptr< definition::SharedElement<E> > DPtr;
    typedef boost::shared_ptr< grid::SharedElement<E> > GPtr;
    typedef std::map<GPtr,DPtr> Map;
    Map unique;
    Grid::ObjectComponentArray::const_iterator gi = input.objects.begin();
    Definition::ObjectComponentSet::iterator di = output.objects.begin();
    for (; gi != input.objects.end(); ++gi, ++di) {
        GPtr gp = gi->template getElement<E>();
        if (!gp) continue;
        std::pair<GPtr,DPtr> item(gp, DPtr());
        std::pair<typename Map::iterator,bool> r = unique.insert(item);
        if (r.second) {
            r.first->second = definition::SharedElement<E>::make(gp->getValue(), gp->isActive());
            if (paramIter != 0 && r.first->second->isActive()) {
                detail::SharedElementTraits<E>::readParameters(
                    paramIter + gp->offset,
                    r.first->second->getValue()
                );
            }
        }
        di->template getElement<E>() = r.first->second;
    }
}

Grid::Grid(Initializer & initializer) : 
    objects(initializer._objects),
    sources(initializer._sources),
    frames(initializer._frames),
    groups(initializer._groups),
    positions(initializer._positions),
    radii(initializer._radii),
    ellipticities(initializer._ellipticities),
    _filterCount(initializer._filterCount),
    _coefficientCount(initializer._coefficientCount),
    _pixelCount(initializer._pixelCount),
    _parameterCount(initializer._parameterCount),
    _constraintCount(initializer._constraintCount),
    _wcs(initializer._wcs)
{
    _filters.swap(initializer._filters);
}

int const Grid::getFilterIndex(FilterId filterId) const {
    FilterMap::const_iterator i = _filters.find(filterId);
    if(i == _filters.end()) {
        throw LSST_EXCEPT(
            lsst::pex::exceptions::InvalidParameterException,
            (boost::format("Filter with ID %d not found.") % filterId).str()
        );
    }
    return i->second;
}

void Grid::writeParameters(ndarray::Array<double,1,1> const & parameters) const {
    double * paramIter = parameters.getData();
    for (PositionArray::const_iterator i = positions.begin(); i != positions.end(); ++i) {
        detail::SharedElementTraits<POSITION>::writeParameters(paramIter + i->offset, i->getValue());
    }
    for (RadiusArray::const_iterator i = radii.begin(); i != radii.end(); ++i) {
        detail::SharedElementTraits<RADIUS>::writeParameters(paramIter + i->offset, i->getValue());
    }
    for (EllipticityArray::const_iterator i = ellipticities.begin(); i != ellipticities.end(); ++i) {
        detail::SharedElementTraits<ELLIPTICITY>::writeParameters(paramIter + i->offset, i->getValue());
    }
}

bool Grid::checkBounds(ndarray::Array<double const,1,1> const & parameters) const {
    double const * paramIter = parameters.getData();
    for (PositionArray::const_iterator i = positions.begin(); i != positions.end(); ++i) {
        if (!i->checkBounds(paramIter)) return false;
    }
    for (RadiusArray::const_iterator i = radii.begin(); i != radii.end(); ++i) {
        if (!i->checkBounds(paramIter)) return false;
    }
    for (EllipticityArray::const_iterator i = ellipticities.begin(); i != ellipticities.end(); ++i) {
        if (!i->checkBounds(paramIter)) return false;
    }
    return true;
}

double Grid::clipToBounds(ndarray::Array<double,1,1> const & parameters) const {
    double * paramIter = parameters.getData();
    double value = 0.0;
    for (PositionArray::const_iterator i = positions.begin(); i != positions.end(); ++i) {
        value += i->clipToBounds(paramIter);
    }
    for (RadiusArray::const_iterator i = radii.begin(); i != radii.end(); ++i) {
        value += i->clipToBounds(paramIter);
    }
    for (EllipticityArray::const_iterator i = ellipticities.begin(); i != ellipticities.end(); ++i) {
        value += i->clipToBounds(paramIter);
    }
    return value;
}

Definition Grid::makeDefinition() const {
    return Initializer::makeDefinition(*this, 0);
}

Definition Grid::makeDefinition(ndarray::Array<double const,1,1> const & parameters) const {
    return Initializer::makeDefinition(*this, parameters.getData());
}

Grid::Ptr Grid::make(Definition const & def) {
    Initializer initializer(def);
    return Grid::Ptr(new Grid(initializer));
}

}}}} // namespace lsst::meas::multifit::grid
