#include "lsst/meas/multifit/grid/Grid.h"
#include <boost/format.hpp>
#include <limits>
#include "lsst/ndarray/eigen.h"

namespace lsst { namespace meas { namespace multifit {

namespace {

template <typename T>
struct DestroyGridElement {
    void operator()(T & p) const { p.~T(); }
};

} // anonymous

namespace grid {

class Initializer {
public:

    template <ParameterType E>
    static void transferElements(
        Grid const & input, Definition & output, double const * paramIter
    ) {
        typedef boost::shared_ptr< definition::ParameterElement<E> > DPtr;
        typedef boost::shared_ptr< grid::ParameterElement<E> > GPtr;
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
                r.first->second = definition::ParameterElement<E>::make(gp->getValue(), gp->isActive());
                if (paramIter != 0 && r.first->second->isActive()) {
                    detail::ParameterElementTraits<E>::readParameters(
                        paramIter + gp->offset,
                        r.first->second->getValue()
                    );
                }
            }
            di->template getElement<E>() = r.first->second;
        }
    }

    template <ParameterType E, typename Container>
    static void transferElements(
        Definition const & input, Grid & output, Container & container, int & offset
    ) {
        typedef boost::shared_ptr< definition::ParameterElement<E> > DPtr;
        typedef boost::shared_ptr< grid::ParameterElement<E> > GPtr;
        typedef std::map<DPtr,GPtr> Map;
        Map unique;
        Definition::ObjectComponentSet::const_iterator di = input.objects.begin();
        Grid::ObjectComponentArray::const_iterator gi = output.objects.begin();
        for (; gi != output.objects.end(); ++gi, ++di) {
            DPtr dp = di->template getElement<E>();
            if (!dp) continue;
            std::pair<DPtr,GPtr> item(dp, GPtr());
            std::pair<typename Map::iterator,bool> r = unique.insert(item);
            if (r.second) {
                if (dp->isActive()) {
                    r.first->second.reset(new grid::ParameterElement<E>(*dp, offset));
                    offset += grid::ParameterElement<E>::SIZE;
                    container._ptrVec.push_back(r.first->second);
                } else {
                    r.first->second.reset(new grid::ParameterElement<E>(*dp, -1));
                }
            }
            GPtr const * gp;
            gi->getElementImpl(gp);
            const_cast<GPtr &>(*gp) = r.first->second;
        }
    }

    static Definition makeDefinition(Grid const & grid, double const * paramIter) {
        Wcs::Ptr wcs;
        if (grid.getWcs()) wcs = grid.getWcs()->clone();
        Definition result(wcs);
        for (Grid::FrameArray::const_iterator i = grid.frames.begin(); i != grid.frames.end(); ++i) {
            result.frames.insert(definition::Frame(*i));
        }
        for (
            Grid::ObjectComponentArray::const_iterator i = grid.objects.begin(); 
            i != grid.objects.end();
            ++i
        ) {
            result.objects.insert(definition::ObjectComponent(*i));
        }
        transferElements<POSITION>(grid, result, paramIter);
        transferElements<RADIUS>(grid, result, paramIter);
        transferElements<ELLIPTICITY>(grid, result, paramIter);
        return result;
    }

    static void initializeGrid(Definition const & input, Grid & output) {
        output.objects._last = output.objects._first 
            = reinterpret_cast<grid::ObjectComponent*>(output._objectData.get());
        output.frames._last = output.frames._first
            = reinterpret_cast<grid::Frame*>(output._frameData.get());
        output.sources._last = output.sources._first
            = reinterpret_cast<grid::SourceComponent*>(output._sourceData.get());
        try {
            int frameCount = 0;
            for (
                Definition::FrameSet::const_iterator i = input.frames.begin();
                i != input.frames.end();
                ++i
            ) {
                std::pair<Grid::FilterMap::iterator,bool> r = output._filters.insert(
                    std::make_pair(i->getFilterId(), output._filterCount)
                );
                if (r.second) ++output._filterCount;
                grid::Frame * newFrame = new (output.frames._last) grid::Frame(
                    *i, output._pixelCount, r.first->second, frameCount
                );
                ++output.frames._last;
                output._pixelCount += newFrame->getPixelCount();
                ++frameCount;
            }
            for (
                Definition::ObjectComponentSet::const_iterator i = input.objects.begin(); 
                i != input.objects.end();
                ++i
            ) {
                grid::ObjectComponent * newObjectComponent = new (output.objects._last) grid::ObjectComponent(
                    *i, output._coefficientCount, frameCount, output._filterCount
                );
                ++output.objects._last;
                output._coefficientCount += newObjectComponent->getCoefficientCount();
                if (newObjectComponent->getBasis()) {
                    _constraintCount += newObjectComponent->getBasis()->getConstraintSize()
                        * (newObjectComponent->getCoefficientCount()
                           / newObjectComponent->getSourceCoefficientCount());
                } else {
                    ++_constraintCount;
                }
            }
            transferElements<POSITION>(input, output, output.positions, output._parameterCount);
            transferElements<RADIUS>(input, output, output.radii, output._parameterCount);
            transferElements<ELLIPTICITY>(input, output, output.ellipticities, output._parameterCount);
            for (
                grid::ObjectComponent * i = output.objects._first; 
                i != output.objects._last;
                ++i
            ) {
                i->validate();
                i->sources._first = output.sources._last;
                for (
                    Grid::FrameArray::const_iterator j = output.frames.begin();
                    j != output.frames.end();
                    ++j
                ) {
                    new (output.sources._last) grid::SourceComponent(*j, *i, output.getWcs());
                    ++output.sources._last;
                }
                i->sources._last = output.sources._last;
            }
        } catch (...) {
            destroyGrid(output);
            throw;
        }
    }

    static void destroyGrid(Grid & g) {
        std::for_each(g.sources._first, g.sources._last, DestroyGridElement<grid::SourceComponent>());
        std::for_each(g.objects._first, g.objects._last, DestroyGridElement<grid::ObjectComponent>());
        std::for_each(g.frames._first, g.frames._last, DestroyGridElement<grid::Frame>());
    }

};

Definition Grid::makeDefinition() const {
    return Initializer::makeDefinition(*this, 0);
}

Definition Grid::makeDefinition(lsst::ndarray::Array<double const,1,1> const & parameters) const {
    return Initializer::makeDefinition(*this, parameters.getData());
}

Grid::Grid(Definition const & def) :
    _filterCount(0),
    _coefficientCount(0),
x    _fluxCoefficientCount(0),
    _morphologyCoefficientCount(0),
    _pixelCount(0),
    _parameterCount(0),
    _constraintCount(0),
    _objectData(new char[sizeof(grid::ObjectComponent) * def.objects.size()]),
    _frameData(new char[sizeof(grid::Frame) * def.frames.size()]),
    _sourceData(new char[sizeof(grid::SourceComponent) * def.frames.size() * def.objects.size()]),
    _wcs()
{
    if (definition.getWcs()) _wcs = definition.getWcs()->clone();
    Initializer::initializeGrid(definition, *this);
}

Grid::~Grid() { Initializer::destroyGrid(*this); }

}}}} // namespace lsst::meas::multifit::grid
